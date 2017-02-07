function [ stimulusStructOut ] = makeImpulseStimStruct( stimulusStruct )
% function [ stimulusStructOut ] = makeImpulseStimStruct( stimulusStruct )
%
% This routine will identify the first time of onset in each column of the
% stimulus struct and then build a stimulus struct that is a set of delta
% functions which are uniformly zero, except for a single value of one at
% the time of onset for each stimulus
%
%   


% set up stimulusStructOut
stimulusStructOut=stimulusStruct;
stimulusStructOut.values(:,:)=0;

% loop through the columns of the stimulus struct and build the impulses
numInstances=size(stimulusStruct.values,1);
for i=1:numInstances

    % Find the stimulus onset so that we can position the impulses. We
    % do that by finding a [0 1] edge from a difference operator.
    tmp = diff(stimulusStruct.values(i,:));
    tmp(tmp < 0) = 0;
    tmp(tmp > 0) = 1;
    
    % Check if the very first value is 1, in which case the stim onset is
    % at the initial value
    if tmp(1)==1
        stimOnset = 1;
    else
        stimOnset = strfind(tmp, [0 1]);
    end
    
    if ~length(stimOnset)==1
        error('Cannot find a unique stimulus onset for this instance')
    end

    stimulusStructOut.values(i,stimOnset)=1;
end

end % function