function [modelResponseStruct] = forwardModelDEDU(params,stimulusStruct)
%% forwardModelDEDU
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of parameters.
%
% The model of neural activity is a step function, controlled by three
% parameters: amplitude, delay (msecs), duration (msecs)
%

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
delayVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'delay'));
durationVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'duration'));


%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(timebaseSecs);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

%% We loop through each column of the stimulus matrix
for i=1:numInstances
    
    % grab the current stimulus
    stimulus=stimulusStruct.values(i,:)';

    % Find the stimulus onsets so that we can adjust the model to it. We
    % do that by finding a [0 1] edge from a difference operator.
    tmp = diff(stimulus');
    tmp(tmp < 0) = 0;
    tmp(tmp > 0) = 1;
    
    % Check if the very first value is 1, in which case the stim onset is
    % at the initial value
    if tmp(1)==1
        stimOnsetIdx=1;
    else
        stimOnsetIdx = strfind(tmp, [0 1]);
    end
    
    % If we can't find a stim onset, return an error.
    if ~length(stimOnsetIdx)==1
        error('Cannot find a unique stimulus onset for this instance')
    end
    
    % Determine the parameters of the step
    stepOnsetTime=stimulusStruct.timebase(stimOnsetIdx)+delayVec(i);
    tempTimeDiff=abs(stimulusStruct.timebase-stepOnsetTime);
    stepOnsetIdx=find(tempTimeDiff==min(tempTimeDiff));
    stepOnsetIdx=stepOnsetIdx(1);
    
    stepOffsetTime=stimulusStruct.timebase(stimOnsetIdx)+delayVec(i)+durationVec(i);
    tempTimeDiff=abs(stimulusStruct.timebase-stepOffsetTime);
    stepOffsetIdx=find(tempTimeDiff==min(tempTimeDiff));
    stepOffsetIdx=stepOffsetIdx(1);
    
    % Make the step
    yNeural=stimulusStruct.values*0;
    yNeural(stepOnsetIdx:stepOffsetIdx)=amplitudeVec(i);
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(i,:)=yNeural;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end