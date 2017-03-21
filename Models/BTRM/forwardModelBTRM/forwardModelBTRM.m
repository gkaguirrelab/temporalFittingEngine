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
%      amplitude - multiplicative scaling of the stimulus.
%      tauGammaIRF - time constant of the neural gamma IRF in msecs. A
%        value of 50 - 100 msecs was found in early visual areas.
%      epsilonCompression - compressive non-linearity of response.
%        Reasonable bouds are [0.1:1]. Not used if dCTS model evoked.
%      tauExpTimeConstant - time constant of the low-pass (exponential
%        decay) component (in mecs). Reasonable bounds [100:100000]
%      nCompression - compressive non-linearity parameter. Reasonable
%        bounds [1:3], where 1 is no compression.
%      divisiveSigma - Adjustment factor to the divisive temporal
%        normalization. Found to be ~0.1 in V1. Set to unity to remove its effect. 
%      tauInhibitoryTimeConstant - time constant of the leaky (exponential)
%        integration of neural signals that produces delayed adaptation.
%      kappaInhibitionAmplitude - multiplicative scaling of the inhibitory
%        effect.

% Hard coded params
vExponent = 3; % Observed by Zaidi et al. to provide the best fit to the RGC response data
dCTS_flag = false;

amplitude_CTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude_CTS'));
tauGammaIRF_CTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauGammaIRF_CTS'));
epsilonCompression_CTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'epsilonCompression_CTS'));
tauInhibitoryTimeConstant_LEAKVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauInhibitoryTimeConstant_LEAK'));
kappaInhibitionAmplitude_LEAKVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'kappaInhibitionAmplitude_LEAK'));
if dCTS_flag
    tauExpTimeConstant_dCTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauExpTimeConstant_dCTS'));
    nCompression_dCTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'nCompression_dCTS'));
    divisiveSigma_dCTSVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'divisiveSigma_dCTS'));
end


%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

% pre-allocate the convolution kernels for speed
gammaKernelStruct.timebase=stimulusStruct.timebase;
gammaKernelStruct.values=stimulusStruct.timebase*0;
inhibitoryKernelStruct.timebase=stimulusStruct.timebase;
inhibitoryKernelStruct.values=stimulusStruct.timebase*0;

if dCTS_flag
    exponentialKernelStruct.timebase=stimulusStruct.timebase;
    exponentialKernelStruct.values=stimulusStruct.timebase*0;
end

%% We loop through each column of the stimulus matrix
for ii=1:numInstances
    
    %% grab the current stimulus
    signalStruct=stimulusStruct;
    signalStruct.values=signalStruct.values(ii,:);
    
    %% Apply gamma convolution
    % Define a gamma function that transforms the
    % stimulus input into a profile of neural activity (e.g., LFP)
    gammaKernelStruct.values = gammaKernelStruct.timebase .* exp(-gammaKernelStruct.timebase/tauGammaIRF_CTSVec(ii));
    % scale the kernel to preserve area of response after convolution
    gammaKernelStruct=normalizeKernelArea(gammaKernelStruct);
    % Convolve the stimulus struct by the gammaKernel
    yNeural=obj.applyKernel(signalStruct,gammaKernelStruct);
    
    %% Implement the dCTS model.
    if dCTS_flag
        % Create the exponential low-pass kernel that defines the time-domain
        % properties of the normalization
        exponentialKernelStruct.values=exp(-1/tauExpTimeConstant_dCTSVec(ii)*exponentialKernelStruct.timebase);
        % scale the kernel to preserve area of response after convolution
        exponentialKernelStruct=normalizeKernelArea(exponentialKernelStruct);
        % Convolve the linear response by the exponential decay
        denominatorStruct=obj.applyKernel(yNeural,exponentialKernelStruct);
        % Apply the compresion and add the semi-saturation constant
        denominatorStruct.values=(divisiveSigma_dCTSVec(ii)^nCompression_dCTSVec(ii)) + ...
            denominatorStruct.values.^nCompression_dCTSVec(ii);
        % Apply the compresion to the numerator
        numeratorStruct.values=yNeural.values.^nCompression_dCTSVec(ii);
        % Compute the final dCTS values
        yNeural.values=numeratorStruct.values./denominatorStruct.values;
    else
        % If we are not implementing the dCTS model, apply compressive
        % non-linearity
        yNeural.values=yNeural.values.^epsilonCompression_CTSVec;
    end
    
    %% Apply amplitude gain
    yNeural.values = yNeural.values.*amplitude_CTSVec(ii);
        
    %% Implement the subtractive effect of a leaky integrator
    % Create an exponential integrator, normalize the area, and convolve
    inhibitoryKernelStruct.values=exp(-1/(tauInhibitoryTimeConstant_LEAKVec(ii)*1000)*inhibitoryKernelStruct.timebase);
    inhibitoryKernelStruct=normalizeKernelArea(inhibitoryKernelStruct);
    inhibitionStruct=obj.applyKernel(yNeural,inhibitoryKernelStruct);
    % raise the inhition temporal profile to vExponent, but preserve area
    inhibitionArea=sum(abs(inhibitionStruct.values));
    inhibitionStruct.values=inhibitionStruct.values.^vExponent;
    inhibitionStruct.values=inhibitionStruct.values/sum(abs(inhibitionStruct.values))*inhibitionArea;
    % Apply the inhibition, scaled by the kappa value
    yNeural.values=yNeural.values- (inhibitionStruct.values*kappaInhibitionAmplitude_LEAKVec(ii));
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(ii,:)=yNeural.values;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end