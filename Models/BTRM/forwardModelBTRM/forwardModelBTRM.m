function [modelResponseStruct] = forwardModelBTRM(obj,params,stimulusStruct)
%% forwardModelBTRM
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of parameters.
%
% The model includes the following stages:
%
%  - a neural impulse response function (modeled as a gamma function)
%  - a compressive non-linearity
%  - a delayed, divisive normalization stage
%    [or, a simple multiplicative exponential decay temporal scaling]
%  - an after-response
%
% The approach is inspired by:
%
%   Zhou, Benson, Kay, Winawer (2016) VSS annual meeting
%   Temporal Summation and Adaptation in Human Visual Cortex
%



%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tauExponentialDecay - time constant of the low-pass (exponential
%        decay) component. Reasonable bounds [0.0001:0.1]
%    These parameters operate at neural timescales, so may be held fixed in
%    the modeling of fMRI data:
%      tauNeuralIRF - time constant of the neural IRF (in msecs). A
%        typical valye might be 100 msecs
%      epsilonCompression - compressive non-linearity parameter. Reasonable
%        bounds [0.1:1]

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
tauExponentialDecayVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauExponentialDecay'));
tauNeuralIRFVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauNeuralIRF'));
epsilonCompressionVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'epsilonCompression'));

% Hard coded params
sigmaValue=0.1;

%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

% pre-allocate the convoluation kernels for speed
gammaKernelStruct.timebase=stimulusStruct.timebase;
gammaKernelStruct.values=stimulusStruct.timebase*0;

exponentialKernelStruct.timebase=stimulusStruct.timebase;
exponentialKernelStruct.values=stimulusStruct.timebase*0;

%% We loop through each column of the stimulus matrix
for ii=1:numInstances
    
    %% grab the current stimulus
    numeratorStruct=stimulusStruct;
    numeratorStruct.values=numeratorStruct.values(ii,:);
            
    %% Apply gamma convolution
    % Define a gamma function that transforms the
    % stimulus input into a profile of neural activity (e.g., LFP)
    gammaKernelStruct.values = gammaKernelStruct.timebase .* exp(-gammaKernelStruct.timebase/tauNeuralIRFVec(ii));
    % scale the kernel to preserve area of response after convolution
    gammaKernelStruct=normalizeKernelAmplitude(gammaKernelStruct);
    % Convolve the stimulus struct by the gammaKernel
    numeratorStruct=obj.applyKernel(numeratorStruct,gammaKernelStruct);
    
    %% Apply amplitude gain
    % scaled by the main response amplitude parameter
    numeratorStruct.values = numeratorStruct.values.*amplitudeVec(ii);

    %% Implement the compressive non-linearity stage
    % Create the exponential low-pass kernel that defines the time-domain
    % properties of the normalization
    exponentialKernelStruct.values=exp(-1*tauExponentialDecayVec(ii)*exponentialKernelStruct.timebase);
    % scale the kernel to preserve area of response after convolution
    exponentialKernelStruct=normalizeKernelAmplitude(exponentialKernelStruct);
    % Convolve the linear response by the exponential decay
    denominatorStruct=obj.applyKernel(numeratorStruct,gammaKernelStruct);
    % Apply the compresion and add the semi-saturation constant
    denominatorStruct.values=(sigmaValue^epsilonCompressionVec(ii)) + ...
        denominatorStruct.values.^epsilonCompressionVec(ii);
    % Apply the compresion to the numerator
    numeratorStruct.values=numeratorStruct.values.^epsilonCompressionVec(ii);
    % Compute the final dCTS values
    yNeural=numeratorStruct.values./denominatorStruct.values;
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(ii,:)=yNeural;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end