function [f,reconstructedTS] = forwardModel(stimMatrix,t,data,paramStruct)

% function [f,reconstructedTS] = forwardModel(stimMatrix,t,data,paramStruct)
%
% takes stimuli and data with a set of parameters, and sees how well those
% parameters allow the model to fit the data

ampVec = paramStruct.Amplitude;
% tau2vec = paramStruct.tau2;

% scale each neural vector by the amplitude parameter, then sum
neuralVec = sum(createNeuralTemporalModelFromStimMatrix(t,stimMatrix,ampVec,paramStruct));

% neural to BOLD
reconstructedTS = neuralVec2BOLD(neuralVec,t,paramStruct.HRF,paramStruct.HRFtimeSamples);

% get error
f = mean((data-reconstructedTS).^2);

gribble = 1;