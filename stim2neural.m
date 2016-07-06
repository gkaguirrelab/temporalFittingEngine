function neuralVec = stim2neural(stimMatrix,t,ampVec,paramStructFixed)

% function neuralVec = stim2neural(stimMatrix,paramMatFit,paramStructFixed)
%
% converts stimulus to neural model

% neural amplitudes
neuralMatrix = repmat(ampVec,[1 size(stimMatrix,2)]).*stimMatrix;

% sum
neuralVec = sum(neuralMatrix);

gribble = 1;