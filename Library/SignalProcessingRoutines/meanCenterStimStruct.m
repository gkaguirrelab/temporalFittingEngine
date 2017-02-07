function [ stimulusStructOut ] = meanCenterStimStruct( stimulusStruct )
% function [ stimulusStructOut ] = meanCenterStimStruct( stimulusStruct )
%
% This routine will mean center each columnn of the stimStruct values
%
%   


% set up stimulusStructOut
stimulusStructOut=stimulusStruct;

% loop through the columns of the stimulus struct and build the impulses
numInstances=size(stimulusStruct.values,1);
for i=1:numInstances

    stimulusStructOut.values(i,:)=stimulusStructOut.values(i,:)-nanmean(stimulusStructOut.values(i,:));
end

end % function