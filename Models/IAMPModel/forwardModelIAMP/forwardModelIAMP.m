function [modelResponseStruct] = forwardModelIAMP(params,stimulusStruct)
%% forwardModelIAMP
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of amplitudes.
%
% This model is extremely simple, and simply scales the stimulus by the
% amplitude.
%

%% Obtain the params
amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));

%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

% We loop through each column of the stimulus matrix
for ii=1:numInstances

    % grab the current stimulus
    stimulus=stimulusStruct.values(ii,:)';
    
    % The neural response is the stimulus input
    % scaled by the amplitude parameter
    responseMatrix(ii,:) = stimulus*amplitudeVec(ii);
        
end % loop over columns of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end % function