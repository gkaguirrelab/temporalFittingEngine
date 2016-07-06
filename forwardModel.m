function [f,reconstructedTS] = forwardModel(stimMatrix,t,data,paramStruct)

% function [f,reconstructedTS] = forwardModel(stimMatrix,t,data,paramStruct)
%
% takes stimuli and data with a set of parameters, and sees how well those
% parameters allow the model to fit the data

paramVecFit = paramStruct.neuralParams;

% scale each neural vector by the amplitude parameter, then sum
neuralVec = stim2neural(stimMatrix,t,paramVecFit,paramStruct);

% neural to BOLD
reconstructedTS = neuralVec2BOLD(neuralVec,t,paramStruct.HRF,paramStruct.HRFtimeSamples);

% get error
f = mean((data-reconstructedTS).^2);

gribble = 1;