function [modelResponseStruct] = forwardModelLEAK(obj, params, stimulusStruct)
%% forwardModelLEAK
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of parameters.
%
% The model includes the following stages:
%
%

%% Hardcoded features

vExponent = 3; % Observed by Zaidi et al. to provide the best fit to the RGC response data
expFilter.timebase = 0:1:30000; % The temporal support of the leaky integrator

%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tau - time constant of leaky exponential integrator that implements
%        the adaptation effect. Reasonable bounds [0.0001:0.1]
%      kappa - multiplicative scaling of the adaptation effect

ampVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
tauVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tau'));
kappaVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'kappa'));

%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

%% We loop through each column of the stimulus matrix
for ii=1:numInstances
    
    % grab the current stimulus
    yNeural.values=stimulusStruct.values(ii,:);
    yNeural.timebase=stimulusStruct.timebase;
    
    %% The neural response begins as the stimulus input
    % scaled by the main response amplitude parameter
    yNeural.values = yNeural.values.*ampVec(ii);
        
    %% Generate the exponential filter
    % This will implement the model of a leaky integrator
    expFilter.values=exp(-1*tauVec(ii)*expFilter.timebase);
    expFilter=normalizeKernelAmplitude(expFilter);
    
    leakyIntegration=obj.applyKernel(yNeural,expFilter);
    yNeural.values=yNeural.values - kappaVec(ii).* (leakyIntegration.values.^3);
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(ii,:)=yNeural.values;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end