function f = forwardModelObjectiveFunction(stimMatrix,t,stimRef,data,paramVecFit,paramStructFixed)

% function f = forwardModelObjectiveFunction(stimMatrix,t,data,paramVec)
%
% same as forwardModel, except it returns just the error. Also, it takes in
% a parameter vector, not a struct, to fit. However, note paramStructFixed,
% which contains parameters that are locked down, like the HRF

% one row of parameters for each stimulus
numberOfParamRows = length(paramStructFixed.stimTag);

% reshape the vector into a matrix 
paramMatFit = reshape(paramVecFit,[numberOfParamRows length(paramVecFit)./numberOfParamRows]);

% scale each neural vector by the amplitude parameter, then sum
neuralVec = stim2neural(stimMatrix,t,stimRef,paramMatFit,paramStructFixed);

% neural to BOLD
reconstructedTS = neuralVec2BOLD(neuralVec,t,paramStructFixed.HRF,paramStructFixed.HRFtimeSamples);

% get error
f = mean((data-reconstructedTS).^2);

gribble = 1;