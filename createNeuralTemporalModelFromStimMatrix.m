function [neuralMatrix] = createNeuralTemporalModelFromStimMatrix(t, stimMatrix, ampVec, paramStructFixed)


%% createNeuralTemporalModelFromStimMatrix
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
%   t - a vector of time points, in seconds
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
%
% 05-30-2016 -  gka wrote it
% 06-23-2016 -  gka modified to serve as a function call in fminsearch
%
%

modelLength = length(t);
stimDimension=size(stimMatrix,1);

%% define default parameters

% parameters of the stimulus
param.afterResponseTiming = 10; % Time after stimulus onset at which the offset response occurs.

% parameters of the neural filters
param.MRamplitude = 1;      % multiplicative scaling of the stimulus into the main neural response. Should be unbounded.
param.ARampRelative = 0;%0.25; % multiplicative scaling of the after-response, relative to main. Reasonable bounds [0:1]
param.tau1 = 0.005;         % time constant of the neural IRF (in seconds). In fMRI data modeling, this will be held fixed.
param.epsilon = .35;        % compressive non-linearity parameter. Reasonable bounds [0.1:1]
param.tau2 = 5;             % time constant of the low-pass (exponential decay) component. Reasonable bounds [0.0011:5]
param.rectify = true;       % controls if rectification is performed upon the neural model.

% currently unused parameters Zhou & Winawer neural model
%param.sigma = 1;        % constant scaling factor of the divisive normalization
%param.n = 0.5;             % the order of the delayed normalization


for s=1:stimDimension
    
    yStimulus=stimMatrix(s,:)';
    
%     offset=yStimulus(1);
%     yStimulus=yStimulus-offset;
    
    param.MRamplitude=ampVec(s);
    
    yNeural = yStimulus.*param.MRamplitude;
%    param.tau2 = tau2vec(s);
    
    %% Implement the neural stage
    % Define the neuralIRF, which is a gamma function that transforms the
    % stimulus input into a profile of neural activity (e.g., LFP)
%    neuralIRF = t .* exp(-t/param.tau1);
    
    % scale to unit sum to preserve amplitude of y following convolution
%    neuralIRF=neuralIRF/sum(neuralIRF);
    
    % Obtain first stage, linear model, which is the stimulus vector with a
    % multiplicative amplitude parameter, then convolved by the neural IRF.
    
%    yLinear = conv(yStimulus*param.MRamplitude,neuralIRF);
    
    % Truncate the convolved vector to the input length. Not sure why I have to
    % do this. Probably mis-using the conv function.
%    yLinear = yLinear(1:modelLength);
    
    %% Implement the compressive non-linearity stage
    % Obtain second stage, CTS model, which is the output of the linear stage
    % subjected to a compressive non-linearity. While this is implemented here
    % as a power law function, it is worth noting that very similar functions
    % are produced by implementing this as an instantaneous divisive
    % normalization.
%    yCTS = yLinear.^param.epsilon;
    
    %% As we are modeling MRI data, we are currently skipping several features
    % of the Zhou & Winawer model. We use the decaying exponential only as a
    % multiplicative weight for the neural response.
    %
    % %% Implement the delayed, divisive normalization stage.
    % % (Instantaneous) divisive normalization is described by:
    % %
    % %  yResponse = Input^n / (sigma^n + Input^n)
    % %
    % % The computation is implemented over time by weighting the input by a
    % % temporal decay function:
    % %
    % %  yResponseDN(t) = Input(t)^n / (sigma^n + (Input(t) * decay(t))^n)
    % %
    % % where '*' indicates here convolution and decay is a decaying exponential
    % % that is defined by the parameter tau2
    % %
    % Create the exponential low-pass function that defines the time-domain
    % properties of the normalization
      decayingExponential=(exp(-1*param.tau2*t));
    %
    % scale to have unit sum to preserve amplitude after convolution
    decayingExponential=decayingExponential/sum(decayingExponential);
    %
    % % convolve the neural response by the decaying exponential
    % yCTSDecayed=conv(yCTS,decayingExponential);
    yNeural=conv(yNeural,decayingExponential);
    %
    % % Truncate the convolved vector to the input length. Not sure why I have to
    % % do this. Probably mis-using the conv function.
    % yCTSDecayed=yCTSDecayed(1:modelLength);
    yNeural=yNeural(1:modelLength);
    %
    % % Assemble the delayed, divisive normalization response
    % yNeuralMR = (yCTS.^param.n) ./ ... % the numerator
    %       ( param.sigma^param.n + yCTSDecayed.^param.n );  % the denominator
    
    %% Implement a decaying exponential scaling of the response
    %
    % This is implemented instead of the divisive normalization approach
    % followed in the original Zhou & Winawer model.
    %
    % The tau2 parameter describes the time constant of this multiplicative
    % scaling.
    %
    
%     % Create the exponential low-pass function that defines the time-domain
%     % properties of the normalization
%     decayingExponential=(exp(-1*param.tau2*t));
%     
%     % Shift and scale the exponential so that the peak coincides with the
%     % initial peak response in the yCTS model. This is needed to unconfound the
%     % overall amplitude response parameter from the exponential decay
%     % parameter.
%     
%     initialPeak=find(yCTS==max(yCTS));
%     initialPeak=initialPeak(1);
%     decayingExponential=circshift(decayingExponential,[0,initialPeak]);
%     decayingExponential=decayingExponential/max(decayingExponential);
%     decayingExponential(1:initialPeak)=1;
    
%     % Apply the exponential function as a mulitiplicative scaling of the
%     % response.
%     
%     yNeuralMR = yCTS.*decayingExponential;
    
    %% Create the after-response
    % This is assumed to be a shifted, scaled version of the main response.
    
%     % Scale, invert, shift, and blank out the circular wrapping)
%     yNeuralAR = yNeuralMR * param.ARampRelative * (-1);
%     yNeuralAR = circshift(yNeuralAR,[0,param.afterResponseTiming]);
%     yNeuralAR(1:param.afterResponseTiming)=NaN;
%     
%     yNeural = nansum([yNeuralMR;yNeuralAR]);
    
    %% Rectify yNeural
    
    % yNeural has a negative after-response, reflecting what is observed in LFP
    % measurements. The param.rectify setting controls if this is converted to
    % a rectified measure to reflect total total synaptic activity, which may
    % be the better model for input to the BOLD system.
    
%     if param.rectify
%         yNeural=abs(yNeural);
%     end
    
    neuralMatrix(s,:)=yNeural - mean(yNeural);
    
end % loop over columns of the stimulus matrix

gribble=1;