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
%  - an adaptive after-response, modeled as the subtractive influence
%    of a leaky (exponentialy decaying) integrator.
%
% The primary, positive response is taken from:
%
%   Zhou, Benson, Kay, Winawer (2017) Systematic changes in temporal
%     summation across human visual cortex
%
% The negative, adaptive effect is taken from:
%
%   Zaidi et al (2012) Neural Locus of Color Afterimages. Current Bio.



%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tauGammaIRF - time constant of the neural gamma IRF in msecs.
%      tauExpTimeConstant - time constant of the low-pass (exponential
%        decay) component (in mecs). Reasonable bounds [100:1000]
%      nCompression - compressive non-linearity parameter. Reasonable
%        bounds [1:3], where 1 is no compression.

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
tauGammaIRFVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauGammaIRF'));
tauExpTimeConstantVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauExpTimeConstant'));
divisiveSigmaVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'divisiveSigma'));
nCompressionVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'nCompression'));
tauInhibitoryTimeConstantVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauInhibitoryTimeConstant'));
kappaInhibitionAmplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'kappaInhibitionAmplitude'));

% Hard coded params
vExponent = 3; % Observed by Zaidi et al. to provide the best fit to the RGC response data


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

inhibitoryKernelStruct.timebase=stimulusStruct.timebase;
inhibitoryKernelStruct.values=stimulusStruct.timebase*0;

%% We loop through each column of the stimulus matrix
for ii=1:numInstances
    
    %% grab the current stimulus
    numeratorStruct=stimulusStruct;
    numeratorStruct.values=numeratorStruct.values(ii,:);
    
    %% Apply gamma convolution
    % Define a gamma function that transforms the
    % stimulus input into a profile of neural activity (e.g., LFP)
    gammaKernelStruct.values = gammaKernelStruct.timebase .* exp(-gammaKernelStruct.timebase/tauGammaIRFVec(ii));
    % scale the kernel to preserve area of response after convolution
    gammaKernelStruct=normalizeKernelArea(gammaKernelStruct);
    % Convolve the stimulus struct by the gammaKernel
    numeratorStruct=obj.applyKernel(numeratorStruct,gammaKernelStruct);
    
    %% Apply amplitude gain
    % scaled by the main response amplitude parameter
    numeratorStruct.values = numeratorStruct.values.*amplitudeVec(ii);
    
    %% Implement the compressive non-linearity stage
    % Create the exponential low-pass kernel that defines the time-domain
    % properties of the normalization
    exponentialKernelStruct.values=exp(-1/tauExpTimeConstantVec(ii)*exponentialKernelStruct.timebase);
    % scale the kernel to preserve area of response after convolution
    exponentialKernelStruct=normalizeKernelArea(exponentialKernelStruct);
    % Convolve the linear response by the exponential decay
    denominatorStruct=obj.applyKernel(numeratorStruct,exponentialKernelStruct);
    % Apply the compresion and add the semi-saturation constant
    denominatorStruct.values=(divisiveSigmaVec(ii)^nCompressionVec(ii)) + ...
        denominatorStruct.values.^nCompressionVec(ii);
    % Apply the compresion to the numerator
    numeratorStruct.values=numeratorStruct.values.^nCompressionVec(ii);
    % Compute the final dCTS values
    yNeural=stimulusStruct;
    yNeural.values=numeratorStruct.values./denominatorStruct.values;
    
    %% Implement the subtractive effect of a leaky integrator
    inhibitoryKernelStruct.values=exp(-1/tauInhibitoryTimeConstantVec(ii)*inhibitoryKernelStruct.timebase);
    % scale the kernel to preserve area of response after convolution
    inhibitoryKernelStruct=normalizeKernelArea(inhibitoryKernelStruct);
    % calculate and apply inhibitory component
    inhibitionStruct=obj.applyKernel(yNeural,inhibitoryKernelStruct);
    % raise the inhition temporal profile to vExponent, but preserve area
    inhibitionArea=sum(abs(inhibitionStruct.values));
    inhibitionStruct.values=inhibitionStruct.values.^vExponent;
    inhibitionStruct.values=inhibitionStruct.values/sum(abs(inhibitionStruct.values))*inhibitionArea;
    yNeural.values=yNeural.values-inhibitionStruct.values*kappaInhibitionAmplitudeVec(ii);
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(ii,:)=yNeural.values;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end