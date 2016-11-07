function [ xValFitStructure, averageResponseStruct, modelResponseStruct ] = crossValidateFits( packetCellArray, tfeHandle, varargin )
% function [ xValFitStructure ] = crossValidateFits( packetCellArray, tfeHandle, varargin )
%
% This is a general cross validation function for the temporalFittingEngine.
%
% Inputs:
%   packetCellArray - a cell array of packets. There must be more than one
%     packet. Each packet must have the same stimLabels and the same
%     unique set of stimTypes. The number and order of stimTypes may vary
%     between the packets.
%  tfeHandle - handle to the model fitting object.
%
% Optional arguments:
%   partitionMethod - how to divide the packets into train and test sets
%     'loo' - leave-one-out for test
%     'split' - split-halves train and test
%     'twentyPercent' - train on 80% of the packets, test on 20%
%     'full' - train and test on all partitions of the packet set
%     'bootstrap' - study only training paritions that are bootstrap
%         resamples of the trainPackets cell array
%   maxPartitions - numeric, specifies maximum number of partitions to
%     evaluate (randomly selected from the available partition set).
%   partitionMatrix - numeric, allows the calling routine to define the
%     partition matrix outside of the cross validation routine.
%   aggregateMethod - method to aggregate paramFits across instances
%     within a packet fit, and to aggregate paramFits and fVals across
%     packets within a partition.
%   verbosity -
%     'none' - the defaut. Shh.
%     'full'
%
% Outputs:
%   xValFitStructure -
%     paramNameCell - the names of the parameters
%     paramMainMatrix - the central tendency of the fit params. The matrix
%       has three dimensions: partition x uniqueStimLabels x param
%     uniqueStimLabels - the names of the stimulus types
%     trainfVals - the vector of fVals from the training partitions
%     testfVals -  - the vector of fVals from the test partitions
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('packetCellArray',@iscell);
p.addRequired('tfeHandle',@(x)(~isempty(x)));
p.addParameter('partitionMethod','loo',@ischar);
p.addParameter('maxPartitions',100,@isnumeric);
p.addParameter('partitionMatrix',[],@isnumeric);
p.addParameter('aggregateMethod','mean',@ischar);
p.addParameter('verbosity','none',@ischar);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.parse(packetCellArray, tfeHandle, varargin{:});

partitionMatrix=p.Results.partitionMatrix;

%% Check the input

% Determine the number of packets
nPackets=length(packetCellArray);

% Error if there is only one packet
if nPackets <= 1
    error('Cross validation requires more than one packet');
end

% check that every packet has the same stimulus labels, and
% handle the case of the stimLabels being strings or numbers
if ischar(packetCellArray{1}.stimulus.metaData.stimLabels{1})
    uniqueStimLabels=unique(packetCellArray{1}.stimulus.metaData.stimLabels);
end
if isnumeric(packetCellArray{1}.stimulus.metaData.stimLabels{1})
    uniqueStimLabels=unique(cell2mat(packetCellArray{1}.stimulus.metaData.stimLabels));
end
for pp=2:nPackets
    if ischar(packetCellArray{1}.stimulus.metaData.stimLabels{1})
        if  ~isequal(uniqueStimLabels, unique(packetCellArray{pp}.stimulus.metaData.stimLabels))
            error('Not all of the packets have the same stimLabels');
        end % test for the same stimLabels
    end % we have string stim labels
    if isnumeric(packetCellArray{1}.stimulus.metaData.stimLabels{1})
        if  ~isequal(uniqueStimLabels, unique(cell2mat(packetCellArray{pp}.stimulus.metaData.stimLabels)))
            error('Not all of the packets have the same stimLabels');
        end % test for the same stimLabels
    end % we have numeric stimLabels
end % loop over the packets


% The partition matrix defines how packets are to be divided into train and
% test sets. It is an n x m matrix, where n is the number of partitions
% over which the cross-validation is to be conducted, and m is the number 
% of packets in the packetCellArray.
% In each row, positive numbers identify packets to be used to train, and
% negative numbers used to test. Packets identified by zeros or NaNs will
% not be used for that partition.
%
% If the paritionMatrix is defined, then the paritionMethod flag is ignored

%% Check a passed partitionMatrix
if ~isempty(partitionMatrix)

    % Check that the second dimension of the matrix is equal to nPackets
    if size(partitionMatrix,2)~=nPackets
        error('There must be one column in the partitionMatrix for each packet');
    end
    
    % Check that there are both negative and positive values
    if sum(partitionMatrix<0)==0 || sum(partitionMatrix>0)==0
        error('The partitionMatrix has all positive or negative values.');
    end
    
    % Reduce the number of partitions if requested
    nPartitions=size(partitionMatrix,1);
    if nPartitions > p.Results.maxPartitions
        % randomly re-assort the rows of the partition matrix, to avoid
        % choosing the same sub-set of limited partitions
        ix=randperm(nPartitions);
        partitionMatrix=partitionMatrix(ix,:);
        partitionMatrix=partitionMatrix(1:p.Results.maxPartitions,:);
        nPartitions=size(partitionMatrix,1);
    end
    
     %% No paritionMatrix was passed, so create it

else
    
    % find all k=2 member partitions of the set of packets
    splitPartitionsCellArray=partitions(1:1:nPackets,2);
    nPartitions=length(splitPartitionsCellArray);
    
    % build a partition matrix, where +=train, -=test
    partitionMatrix=zeros(nPartitions*2,nPackets);
    for ss=1:nPartitions
        partitionMatrix(ss,splitPartitionsCellArray{ss}{1})=1;
        partitionMatrix(ss,splitPartitionsCellArray{ss}{2})=-1;
        partitionMatrix(end+1-ss,splitPartitionsCellArray{ss}{2})=1;
        partitionMatrix(end+1-ss,splitPartitionsCellArray{ss}{1})=-1;
    end % loop over the partitions
    
    % restrict the partition Matrix following partitionMethod
    switch p.Results.partitionMethod
        case 'loo'
            nTargetTrain=nPackets-1;
            partitionsToUse=find(sum(partitionMatrix>0,2)==nTargetTrain);
            partitionMatrix=partitionMatrix(partitionsToUse,:);
        case 'splitHalf'
            nTargetTrain=floor(nPackets/2);
            partitionsToUse=find(sum(partitionMatrix>0,2)==nTargetTrain);
            partitionMatrix=partitionMatrix(partitionsToUse,:);
        case 'twentyPercent'
            nTargetTrain=ceil(nPackets*0.8);
            partitionsToUse=find(sum(partitionMatrix>0,2)==nTargetTrain);
            partitionMatrix=partitionMatrix(partitionsToUse,:);
        case 'full'
            partitionMatrix=partitionMatrix;
        case 'bootstrap'
            partitionMatrix=[];
        otherwise
            error('Not a recognized partion Method');
    end % switch on partition method
    
    % Reduce the number of partitions if requested
    nPartitions=size(partitionMatrix,1);
    if nPartitions > p.Results.maxPartitions
        % randomly re-assort the rows of the partition matrix, to avoid
        % choosing the same sub-set of limited partitions
        ix=randperm(nPartitions);
        partitionMatrix=partitionMatrix(ix,:);
        partitionMatrix=partitionMatrix(1:p.Results.maxPartitions,:);
        nPartitions=size(partitionMatrix,1);
    end % check for maximum desired number of partitions
end % Check or create partitionMatrix

%% Loop through the packetCellArray and obtain the fits

% Obtain the uniqueStim types for the first packet. We will make sure this
% is the same across all packets
uniqueStimTypes=unique(packetCellArray{1}.stimulus.metaData.stimTypes);

for pp=1:nPackets
    
    % Check that every train and test packet have the same unique stimTypes
    if pp>1
        if ~isequal(uniqueStimTypes, unique(packetCellArray{pp}.stimulus.metaData.stimTypes))
            error('The packets have different unique stimTypes.');
        end % test for the same stimTypes
    end % not on the first packet
    
    % report our progress
    if strcmp(p.Results.verbosity,'full')
        fprintf('* Fitting packet <strong>%g</strong> of <strong>%g</strong>\n', pp, nPackets);
    end
    
    % Get the number of instances in this packet. If the
    % defaultParamsInfo is defined, use this. If not, check to see if
    % the field stimulus.metaData.nInstances is defined. If not, take
    % the number of rows in stimulus.values as the number of instances.
    if ~isempty(p.Results.defaultParamsInfo)
        defaultParamsInfo.nInstances=p.Results.defaultParamsInfo.nInstances;
    else
        if isfield(packetCellArray{pp}.stimulus.metaData,'nInstances')
            defaultParamsInfo.nInstances=packetCellArray{pp}.stimulus.metaData.nInstances;
        else
            defaultParamsInfo.nInstances=size(packetCellArray{pp}.stimulus.values,1);
        end % check for stimulus.metaData.nInstances
    end % check for defaultParamsInfo
    
    % do the fit
    [paramsFit,fVal,~] = ...
        tfeHandle.fitResponse(packetCellArray{pp},...
        'defaultParamsInfo', defaultParamsInfo, ...
        varargin{:});
    
    % Aggregate parameters across instances by stimulus types
    for cc=1:length(uniqueStimTypes)
        thisStimTypeIndices=find(packetCellArray{pp}.stimulus.metaData.stimTypes==uniqueStimTypes(cc));
        switch p.Results.aggregateMethod
            case 'mean'
                eachPacketParamsFit(pp,cc,:)=mean(paramsFit.paramMainMatrix(thisStimTypeIndices,:),1);
            case 'median'
                eachPacketParamsFit(pp,cc,:)=median(paramsFit.paramMainMatrix(thisStimTypeIndices,:),1);
            otherwise
                error('This is an undefined aggregation method');
        end % switch over aggregation methods
    end % loop over uniqueStimLabels
    
    % save the fVal for this trainPacket fit
    eachPacketfVals(pp)=fVal;
end


%% Loop through the partitions
testfVals=[];
for pp=1:nPartitions
    
    if strcmp(p.Results.partitionMethod,'bootstrap')
        trainIdx=randsample(1:nPackets,nPackets,true);
        testIdx=[];
    else
        partitionRow=partitionMatrix(pp,:);
        trainIdx=find(partitionRow>0);
        testIdx=find(partitionRow<0);
    end
        
    % Empty the variables that aggregate across the trainPackets
    
    % Aggregate across trainining packets
    switch p.Results.aggregateMethod
        case 'mean'
            % detect the case of a single parameter type
            if ndims(eachPacketParamsFit)==2
                trainParamsFit(pp,:,1)=mean(eachPacketParamsFit(trainIdx,:),1);
            end
            if ndims(eachPacketParamsFit)==3
                trainParamsFit(pp,:,:)=mean(eachPacketParamsFit(trainIdx,:,:),1);
            end
            trainfVals(pp)=mean(eachPacketfVals(trainIdx));
        case 'median'
            if ndims(eachPacketParamsFit)==2
                trainParamsFit(pp,:,1)=median(eachPacketParamsFit(trainIdx,:),1);
            end
            if ndims(eachPacketParamsFit)==3
                trainParamsFit(pp,:,:)=median(eachPacketParamsFit(trainIdx,:,:),1);
            end
            trainfVals(pp)=median(eachPacketfVals(trainIdx));
        otherwise
            error('This is an undefined aggregation method');
    end % switch over aggregation methods
    
    % Check if there are testPackets to test
    if ~isempty(testIdx)
        testfValsLocal=[];

        % Loop over the test set
        for tt=1:length(testIdx)
            
            % report our progress
            if strcmp(p.Results.verbosity,'full')
                fprintf('* Partition, test packet <strong>%g</strong> , <strong>%g</strong>\n', pp, tt);
            end

            % get the number of instances in this packet
            if isempty(p.Results.defaultParamsInfo)
                defaultParamsInfo.nInstances=size(packetCellArray{testIdx(tt)}.stimulus.values,1);
            else
                defaultParamsInfo.nInstances=p.Results.defaultParamsInfo.nInstances;
            end
            
            % Build a params structure that holds the prediction from the train set
            predictParams = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
            for ii=1:defaultParamsInfo.nInstances
                predictParams.paramMainMatrix(ii,:)= trainParamsFit(pp, packetCellArray{testIdx(tt)}.stimulus.metaData.stimTypes(ii), :);
            end
            
            % create a paramsVec to pass to the fitError method
            predictParamsVec=tfeHandle.paramsToVec(predictParams);
            
            % get the error for the prediction of this test packet
            [fVal,~] = tfeHandle.fitError(predictParamsVec,packetCellArray{testIdx(tt)},varargin{:});
            
            % save the fVal for this trainPacket fit
            testfValsLocal(tt)=fVal;
            
        end % loop over instances in this test packet
        
        % Aggregate across test packets
        switch p.Results.aggregateMethod
            case 'mean'
                testfVals(pp)=mean(testfValsLocal);
            case 'median'
                testfVals(pp)=median(testfValsLocal);
            otherwise
                error('This is an undefined aggregation method');
        end % switch over aggregation methods
    end % check that there are testPackets
end % loop over partitions

% Assemble the xValFitStructure for return
tmpParams = tfeHandle.defaultParams();
xValFitStructure.paramNameCell=tmpParams.paramNameCell;
xValFitStructure.uniqueStimLabels=uniqueStimLabels;
xValFitStructure.paramMainMatrix=trainParamsFit;
xValFitStructure.trainfVals=trainfVals;
xValFitStructure.testfVals=testfVals;

% Create an average signal and an average fit, using central tendency of
% the fits from the train packets, and the response.values from all
% packets. This can be performed only if the train and test packets have
% the equal stimTypes, and equal response timebases.
averageResponseStruct=[];
modelResponseStruct=[];
stimTypes=packetCellArray{1}.stimulus.metaData.stimTypes;
respTimebase=packetCellArray{1}.response.timebase;
for pp=1:nPackets
    if ~isempty(packetCellArray{pp})
        checkStimTypes(pp)=isequal(stimTypes,packetCellArray{pp}.stimulus.metaData.stimTypes);
        checkResponseTimebase(pp)=isequal(respTimebase,packetCellArray{pp}.response.timebase);
    end
end
if sum(checkStimTypes)==nPackets && ...
        sum(checkResponseTimebase)==nPackets
    
    % Save the response timebase
    averageResponseStruct.timebase=respTimebase;
    
    % Obtain the central tendency of the response.values
    for pp=1:nPackets
        responseMatrix(pp,:)=packetCellArray{pp}.response.values;
    end % loop over comboPackets
    
    switch p.Results.aggregateMethod
        case 'mean'
            averageResponseStruct.values=mean(responseMatrix);
            averageResponseStruct.metaData.sem=std(responseMatrix)/sqrt(size(responseMatrix,1));
        case 'median'
            averageResponseStruct.values=median(responseMatrix);
        otherwise
            error('This is an undefined aggregation method');
    end % switch over aggregation methods
    
    % Obtain the model fit testPackets were run
    if ~isempty(testfVals)
        
        % get the number of instances in the first packet
        if isempty(p.Results.defaultParamsInfo)
            defaultParamsInfo.nInstances=size(packetCellArray{1}.stimulus.values,1);
        else
            defaultParamsInfo.nInstances=p.Results.defaultParamsInfo.nInstances;
        end
        
        % Obtain the central tendency of the parameters found in the
        % training set
        switch p.Results.aggregateMethod
            case 'mean'
                averageTrainParams=mean(trainParamsFit,1);
            case 'median'
                averageTrainParams=median(trainParamsFit,1);
            otherwise
                error('This is an undefined aggregation method');
        end % switch over aggregation methods
        
        % Build a params structure that holds the prediction from the training set
        predictParams = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
        for ii=1:defaultParamsInfo.nInstances
            predictParams.paramMainMatrix(ii,:)= averageTrainParams(1, stimTypes(ii), :);
        end
        
        % get the error for the prediction of this test packet
        modelResponseStruct = tfeHandle.computeResponse(predictParams,comboPackets{1}.stimulus,comboPackets{1}.kernel,'AddNoise',false);
    end % Check if there are any testPackets
end % all stimTypes are the same, so can build the responseStructs

end % function

