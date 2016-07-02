function neuralVec = stim2neural(stimMatrix,t,stimRef,paramVecFit,paramStructFixed)

% function neuralVec = stim2neural(stimMatrix,paramMatFit,paramStructFixed)
%
% converts stimulus to neural model

% multiply assignment matrix with parameter vector
ampMatrix = stimRef*paramVecFit;

% neural amplitudes
neuralMatrix = repmat(ampMatrix,[1 size(stimMatrix,2)]).*stimMatrix;

% sum
neuralVec = sum(neuralMatrix);

gribble = 1;