function [modelResponseStruct] = forwardModelDEDU(params,stimulusStruct)
%% forwardModelDEDU
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of parameters.
%
% The model of neural activity is a step function, controlled by three
% parameters: amplitude, delay (secs), duration (secs)
%

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
delayVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'delay'));
durationVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'duration'));


%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% the model has parameters that are tuned for units of seconds, so
% we convert our timebase
timebaseMsecs=stimulusStruct.timebase;
timebaseSecs=timebaseMsecs/1000;

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
    
    % If the stimOnsetIdx is empty, check if the first value of the
    % stimulus vector is different from the last. This tells us that the
    % stimulus begins at the first timepoint
    if (isempty(stimOnsetIdx) && ~(stimulus(1)==stimulus(end)))
        stimOnsetIdx=1;
    end
    
    % If we can't find a stim onset, return an error.
    if ~length(stimOnsetIdx)==1
        error('Cannot find a unique stimulus onset for this instance. Assuming start at first time point.')
    end
    
    % Determine the parameters of the step
    stepOnsetTime=timebaseSecs(stimOnsetIdx)+delayVec(i);
    tempTimeDiff=abs(timebaseSecs-stepOnsetTime);
    stepOnsetIdx=find(tempTimeDiff==min(tempTimeDiff));
    stepOnsetIdx=stepOnsetIdx(1);
    
    stepOffsetTime=timebaseSecs(stimOnsetIdx)+delayVec(i)+durationVec(i);
    tempTimeDiff=abs(timebaseSecs-stepOffsetTime);
    stepOffsetIdx=find(tempTimeDiff==min(tempTimeDiff));
    stepOffsetIdx=stepOffsetIdx(1);
    
    % Make the step
    yNeural=timebaseSecs*0;
    yNeural(stepOnsetIdx:stepOffsetIdx)=amplitudeVec(i);
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(i,:)=yNeural;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end