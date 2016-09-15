function [neuralMatrix] = forwardModel(t, stimMatrix, ampVec)
%% createNeuralResponseFromStimMatrix
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of amplitudes.
%
% This model is extremely simple, and simply scales the stimulus by the
% amplitude.
%
%
% Input properties:
%
%   t - a vector of time points, in milliseconds. It is unused in this
%   model.
%   stimMatrix - DESCRIBE THIS
%   ampVec - a vector of amplitudes, one for each instance
%
% Output properties:
%
%   neuralMatrix
%
% 09-13-2016 -  gka wrote it

% determine how many stimulus instances there are
stimDimension=size(stimMatrix,1);

% We loop through each column of the stimulus matrix
for instance=1:stimDimension
        
    %% The neural response is the stimulus input
    % scaled by the amplitude parameter
    neuralMatrix(instance,:) = stimMatrix(instance,:).*ampVec(instance);
        
end % loop over columns of the stimulus matrix

end