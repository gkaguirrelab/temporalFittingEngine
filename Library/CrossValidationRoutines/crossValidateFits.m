function [ xValFitStructure ] = crossValidateFits( packetCellArray, tfeHandle, varargin )
%function [ xValFitStructure ] = crossValidateFits( packetCellArray, tfeHandle, varargin )
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
%   maxPartitions - numeric, specifies maximum number of partitions to
%     evaluate (randomly selected from the available partition set).
%   aggregateMethod - method to aggregate paramFits across instances
%     within a packet fit, and to aggregate paramFits and fVals across
%     packets within a partition.
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
p.addParameter('aggregateMethod','mean',@ischar);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.parse(packetCellArray, tfeHandle, varargin{:});

%% Check the input

% Determine the number of packets
nPackets=length(packetCellArray);

% Error if there is only one packet
if nPackets <= 1
    error('cross validation requires more than one packet');
end

% check that every packet has the same stimulus labels and the same unique
% set of stimulus types
uniqueStimTypes=unique(packetCellArray{1}.stimulus.metaData.stimTypes);
uniqueStimLabels=unique(packetCellArray{1}.stimulus.metaData.stimLabels);
if ~(length(uniqueStimTypes)==length(uniqueStimLabels))
    error('There must be as many unique stimTypes as stimLabels');
end
for pp=2:nPackets
    if ~isequal(uniqueStimTypes, unique(packetCellArray{pp}.stimulus.metaData.stimTypes))
        error('Not all of the packets have the same stimTypes');
    end % test for the same stimTypes
    if  ~isequal(uniqueStimLabels, unique(packetCellArray{pp}.stimulus.metaData.stimLabels))
        error('Not all of the packets have the same stimLabels');
    end % test for the same stimLabels
end

%% Calculate a partition matrix
% find all k=2 member partitions of the set of packets
splitPartitionsCellArray=partitions(1:1:nPackets,2);
nPartitions=length(splitPartitionsCellArray);

% build a partition matrix, where 1=train, 0=test
partitionMatrix=zeros(nPartitions*2,nPackets);
for ss=1:nPartitions
    partitionMatrix(ss,splitPartitionsCellArray{ss}{1})=1;
    partitionMatrix(end+1-ss,splitPartitionsCellArray{ss}{2})=1;
end % loop over the partitions

% restrict the partition Matrix following partitionMethod
switch p.Results.partitionMethod
    case 'loo'
        nTargetTrain=nPackets-1;
        partitionsToUse=find(sum(partitionMatrix,2)==nTargetTrain);
        partitionMatrix=partitionMatrix(partitionsToUse,:);
    case 'splitHalf'
        nTargetTrain=floor(nPackets/2);
        partitionsToUse=find(sum(partitionMatrix,2)==nTargetTrain);
        partitionMatrix=partitionMatrix(partitionsToUse,:);
    case 'twentyPercent'
        nTargetTrain=ceil(nPackets*0.8);
        partitionsToUse=find(sum(partitionMatrix,2)==nTargetTrain);
        partitionMatrix=partitionMatrix(partitionsToUse,:);
    case 'full'
        partitionMatrix=partitionMatrix;
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

%% Loop through the partitions of the packetCellArray and fit
% Empty the variables that aggregate across the partitions
trainParamsFit=[];
trainfVals=[];
testfVals=[];

for pp=1:nPartitions
    trainPackets = packetCellArray(partitionMatrix(pp,:)==1);
    testPackets = packetCellArray(partitionMatrix(pp,:)==0);
    
    % Empty the variables that aggregate across the trainPackets
    trainParamsFitLocal=[];
    trainfValsLocal=[];
    
    % Loop through the trainPackets
    for tt=1:length(trainPackets)
        
        % report our progress
        fprintf('* Partition, packet <strong>%g</strong> , <strong>%g</strong>\n', pp, tt);
        
        % get the number of instances in this packet
        if isempty(p.Results.defaultParamsInfo)
            defaultParamsInfo.nInstances=size(trainPackets{tt}.stimulus.values,1);
        else
            defaultParamsInfo.nInstances=p.Results.defaultParamsInfo.nInstances;
        end
        
        % do the fit
        [paramsFit,fVal,~] = ...
            tfeHandle.fitResponse(trainPackets{tt},...
            'defaultParamsInfo', defaultParamsInfo, ...
            varargin{:});
        
        % Aggregate parameters across instances by stimulus types
        for cc=1:length(uniqueStimTypes)
            thisStimTypeIndices=find(trainPackets{tt}.stimulus.metaData.stimTypes==uniqueStimTypes(cc));
            switch p.Results.aggregateMethod
                case 'mean'
                    trainParamsFitLocal(tt,cc,:)=mean(paramsFit.paramMainMatrix(thisStimTypeIndices,:),1);
                case 'mean'
                    trainParamsFitLocal(tt,cc,:)=median(paramsFit.paramMainMatrix(thisStimTypeIndices,:),1);
                otherwise
                    error('This is an undefined aggregation method');
            end % switch over aggregation methods
        end % loop over uniqueStimLabels
        
        % save the fVal for this trainPacket fit
        trainfValsLocal(tt)=fVal;
        
    end % loop over trainPackets
    
    % Aggregate across trainining packets
    switch p.Results.aggregateMethod
        case 'mean'
            trainParamsFit(pp,:,:)=mean(trainParamsFitLocal,1);
            trainfVals(pp)=mean(trainfValsLocal);
        case 'median'
            trainParamsFit(pp,:,:)=median(trainParamsFitLocal,1);
            trainfVals(pp)=median(trainfValsLocal);
        otherwise
            error('This is an undefined aggregation method');
    end % switch over aggregation methods
    
    % loop over the test set
    for tt=1:length(testPackets)
        
        % get the number of instances in this packet
        if isempty(p.Results.defaultParamsInfo)
            defaultParamsInfo.nInstances=size(testPackets{tt}.stimulus.values,1);
        else
            defaultParamsInfo.nInstances=p.Results.defaultParamsInfo.nInstances;
        end
        
        % Build a params structure that holds the prediction from the test set
        predictParams = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
        for ii=1:defaultParamsInfo.nInstances
            predictParams.paramMainMatrix(ii,:)= trainParamsFit(pp, testPackets{tt}.stimulus.metaData.stimTypes(ii), :);
        end
        
        % create a paramsVec to pass to the fitError method
        predictParamsVec=tfeHandle.paramsToVec(predictParams);
        
        % get the error for the prediction of this test packet
        [fVal,~] = tfeHandle.fitError(predictParamsVec,testPackets{tt},varargin{:});
        
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
    
end % loop over partitions

% Assemble the xValFitStructure for return
xValFitStructure.uniqueStimLabels=uniqueStimLabels;
xValFitStructure.paramNameCell=predictParams.paramNameCell;
xValFitStructure.paramMainMatrix=trainParamsFit;
xValFitStructure.trainfVals=trainfVals;
xValFitStructure.testfVals=testfVals;

end % function

