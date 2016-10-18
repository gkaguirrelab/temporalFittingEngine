function [ trainParamsFit, trainfVals ] = crossValidateFits( packetCellArray, tfeHandle, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('packetCellArray',@iscell);
p.addRequired('tfeHandle',@(x)(strcmp(class(x),'function_handle')));
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('aggregateMethod','mean',@ischar);
p.parse(paramsCellArray, tfeHandle, varargin{:});

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
% Currently a dumb LOO
nPartitions=nPackets;
testSets=eye(nPackets);

% Empty the variables that aggregate across the partitions.
trainParamsFit=[];
trainfVals=[];
testfVals=[];

%% Loop through the partitions of the packetCellArray and fit
for pp=1:nPartitions
    trainPackets = packetCellArray(testSets(pp,:)==0);
    testPackets = packetCellArray(testSets(pp,:)==1);
    
    % Empty the variables that aggregate across the trainPackets
    trainParamsFitLocal=[];
    trainfValsLocal=[];
    
    % Loop through the trainPackets
    for tt=1:length(trainPackets)
        
        % get the number of instances in this packet
        if isempty(defaultParamsInfo)
            defaultParamsInfo.nInstances=size(trainPackets{tt}.stimulus.values,1);
        end
        
        % do the fit
        [paramsFit,fVal,~] = ...
            tfeHandle.fitResponse(thePacket,...
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
            trainParamsFit(pp,:,:)=mean(trainParamsFit,1);
            trainfVals(pp)=mean(trainfValsLocal);
        case 'median'
            trainParamsFit(pp,:,:)=median(trainParamsFit,1);
            trainfVals(pp)=median(trainfValsLocal);
        otherwise
            error('This is an undefined aggregation method');
    end % switch over aggregation methods
    
    % Apply the model to the test data
    
end % loop over partitions


end % function

