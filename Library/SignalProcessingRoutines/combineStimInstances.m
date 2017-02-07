function [ stimulusStructOut ] = combineStimInstances( stimulusStruct )
% function [ stimulusStructOut ] = combineStimInstances( stimulusStruct )
%
% This routine will identify the those stimulus instances that have the
% same stimType and add them together.
%

% Check that the metaData contains stimTypes and stimLabels
if ~isfield(stimulusStruct.metaData,'stimLabels') || ~isfield(stimulusStruct.metaData,'stimTypes')
    error('The stimulusStruct must have stimLabels and stimTypes defined');
end % check the first packet for stimTypes

uniqueStimTypes=unique(stimulusStruct.metaData.stimTypes);

% set up stimulusStructOut
stimulusStructOut=stimulusStruct;
stimulusStructOut.values=zeros(length(uniqueStimTypes),length(stimulusStruct.timebase));

% loop through the unique stimTypes and combine the matching stimulus
% instances
for i=1:length(uniqueStimTypes)
    idx=find(stimulusStruct.metaData.stimTypes==uniqueStimTypes(i));
    stimulusStructOut.values(i,:)=nansum(stimulusStruct.values(idx,:));
end

stimulusStructOut.metaData.stimTypes=uniqueStimTypes;

end % function