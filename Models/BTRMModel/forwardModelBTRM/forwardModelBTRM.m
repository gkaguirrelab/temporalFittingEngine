function [modelResponseStruct] = forwardModelBTRM(params,stimulusStruct)
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
%   With additional modeling details taken from:
%
%   McLelland, D., Ahmed, B., & Bair, W. (2009). Responses to static visual
%   images in macaque lateral geniculate nucleus: implications for
%   adaptation, negative afterimages, and visual fading.
%   The Journal of Neuroscience, 29(28), 8996-9001.
%
%   McLelland, D., Baker, P. M., Ahmed, B., & Bair, W. (2010).
%   Neuronal responses during and after the presentation of static visual
%   stimuli in macaque primary visual cortex.
%   The Journal of Neuroscience, 30(38), 12619-12631.
%
% This model implements a negative after-response. Based upon the findings
%   of McLelland et al., 2010, it seems that (for the most part) V1 neurons
%   respond with an after-response that has an amplitude that is
%   proportional to the main response, with a very similar time constant.
%
%   We adopt these properties in the model by treating the after-response
%   as an inverted, shifted, and scaled version of the main-response.
%
% Input properties:
%
%   t - a vector of time points, in milliseconds
%   yStimulus - a vector of stimulus amplitude. The length must be the same
%               as t. The absolute amplitude is arbitraty.
%   displayFitPlotIn - Boolean flag indicating if you want a plot. Optional.
%   paramIn - a structure of parameters that define the model. These are
%             described below. Optional.
%
% Output properties:
%
%   yBOLD - a vector of response amplitudes, of the same length as t.
%   t - the temporal vector
%   paramOut - the parameters passed back out
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

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
tauExponentialDecayVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauExponentialDecay'));
tauNeuralIRFVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'tauNeuralIRF'));
epsilonCompressionVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'epsilonCompression'));

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

%% We loop through each column of the stimulus matrix
for i=1:numInstances
    
    % grab the current stimulus
    stimulus=stimulusStruct.values(i,:)';
    
    %% The neural response begins as the stimulus input
    % scaled by the main response amplitude parameter
    yNeural = stimulus*amplitudeVec(i);
    
    %% find the initial peak of the scaled stimulus
    initialPeakPoint=find(abs(yNeural)==max(abs(yNeural)));
    initialPeakPoint=initialPeakPoint(1);
    initialPeakValue=yNeural(initialPeakPoint);
    
    %% Apply gamma convolution
    % Define a gamma function that transforms the
    % stimulus input into a profile of neural activity (e.g., LFP)
    gammaIRF = timebaseSecs .* exp(-timebaseSecs/tauNeuralIRFVec(i));
    
    %scale the IRF to preserve area of response after convolution
    gammaIRF=gammaIRF./sum(gammaIRF);
    
    % Obtain first stage, linear model, which is the scaled stimulus
    % convolved by the neural IRF.
    yNeural = conv(yNeural,gammaIRF);
    yNeural = yNeural(1:modelLength);
    
    %% Implement the compressive non-linearity stage
    % Obtain second stage, CTS model, which is the output of the linear stage
    % subjected to a compressive non-linearity. While this is implemented here
    % as a power law function, it is worth noting that very similar functions
    % are produced by implementing this as an instantaneous divisive
    % normalization.
    yNeural = yNeural.^epsilonCompressionVec(i);
    
    % Restore the peak signed, abs amplitude
    yNeural=(yNeural/yNeural(initialPeakPoint))*initialPeakValue;
    
    % Create the exponential low-pass function that defines the time-domain
    % properties of the normalization
    decayingExponential=(exp(-1*tauExponentialDecayVec(i)*timebaseSecs));
    
    % Position the decaying exponential to have unit value at that time
    % point of the initial peak.
    decayingExponential=circshift(decayingExponential,[0,initialPeakPoint]);
    decayingExponential=decayingExponential/max(decayingExponential);
    decayingExponential(1:initialPeakPoint)=1;
    
    % Apply the exponential decay as a multiplicative scaling
    yNeural=yNeural.*decayingExponential;
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(s,:)=yNeural;
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end