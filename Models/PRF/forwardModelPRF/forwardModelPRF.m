function [modelResponseStruct] = forwardModelPRF(params,stimulusStruct)
%% forwardModelPRF
%
% This function creates a model of neural response given a movie of
% stimulus input, a vector of time points, and a set of parameters.
%
%
% The approach is inspired by:
%


%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tauExponentialDecay - time constant of the low-pass (exponential
%        decay) component. Reasonable bounds [0.0001:0.1]
%    These parameters operate at neural timescales, so may be held fixed in
%    the modeling of fMRI data:
%      tauNeuralIRF - time constant of the neural IRF (in seconds). A
%        typical valye might be 0.005 secs (5 msecs)
%      epsilonCompression - compressive non-linearity parameter. Reasonable
%        bounds [0.1:1]

xPositionVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'xPosition'));
yPositionVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'yPosition'));
sigmaSizeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sigmaSize'));

%% Define basic model features

% the model has parameters that are tuned for units of seconds, so
% we convert our timebase
timebaseMsecs=stimulusStruct.timebase;
timebaseSecs=timebaseMsecs/1000;

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(timebaseSecs);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

%% Perform the forward model calculation

% MAKE THIS CODE FIT THE VARIABLE NAMES HERE

% Gaussian windowed sum of the stimulus at the x and y position.
[meshX , meshY] = meshgrid(1:size(stim,2),1:size(stim,1));
f = exp (-((meshY-y0).^2 + (meshX-x0).^2) ./ (2*s.^2));
for i =1:size(stim,3)
    S = stim(:,:,i).*f;
    gStim(i) = sum(S(:));
end

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end