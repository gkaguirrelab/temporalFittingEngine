
% NEED TO IMPLEMENT AN ADJUSTMENT OF PARAMETERS IN EQUATIONS TO ACCOUNT FOR
% SAMPLING RESOLUTION OF MODEL

close all;

% parameters of the neural filters

param.tau1 = 5/1000; % time constant of the neural IRF (in seconds)
param.epsilon = 0.35; % compressive non-linearity parameter
param.tau2 = 0.5; % time constant of the low-pass component of delayed normalization
param.sigma = 1; % parameter of the delayed normalization
param.n = 2; % the order of the delayed normalization

% parameters of the hemodynamic filters

param.gamma1 = 6;
param.gamma2 = 9;
param.surround = 4;

% parameters of the stimulus simulation

modelDuration=40; % total duration of the modeled period, in seconds
modelResolution=1; % temporal sampling frequency of model, in Hz
modelLength=modelDuration*modelResolution; % duh
stimulusDuration=12; % duration of the modeled step of neural activity, in seconds

% Create a stimulus model, which is a simple step function of input

t = linspace(1,modelDuration,modelDuration*modelResolution); % the temporal domain
yStimulus = t*0; % the stimulus vector
yStimulus(2:(stimulusDuration*modelResolution)+1) = 1; % implement the step function of input

% Define the neuralIRF, which is a gamma function that transforms the
% stimulus input into a profile of neural activity (e.g., LFP)

neuralIRF = t .* exp(-t/param.tau1); 
neuralIRF=neuralIRF/sum(neuralIRF); % scale to unit sum to preserve
                                    % amplitude of y following convolution

% Obtain first stage, linear model, which is the stimulus vector convolved
% by the neural IRF

yLinear = conv(yStimulus,neuralIRF);
yLinear = yLinear(1:modelLength); % Need to explore why I have to truncate
                                  % what is returned by the conv function
                                  
% Obtain second stage, CTS model, which is the output of the linear stage
% subjected to a compressive non-linearity

yCTS = yLinear.^param.epsilon;

% Implement the delayed normalization

% Create the exponential low-pass function that defines the time-domain
% properties of the normalization

decayingExponential=(exp(-1*param.tau2*t));
decayingExponential=decayingExponential/sum(decayingExponential);
lowpassNeuralFilter=conv(yCTS,decayingExponential);
lowpassNeuralFilter=lowpassNeuralFilter(1:modelLength); % Need to explore why I have to truncate
                                  % what is returned by the conv function
                                  
% Assemble the delayed normalization response

yDN = (yCTS.^(param.n)) ./ ... % the numerator 
      ( param.sigma^param.n + lowpassNeuralFilter.^param.n );

% Create a double-gamma model of the hemodynamic response function (HRF)

BOLDHRF = gampdf(t, param.gamma1, 1) - gampdf(t, param.gamma2, 1)/param.surround;
BOLDHRF = BOLDHRF/sum(BOLDHRF);  % scale to unit sum to preserve
                                 % amplitude of y following convolution

% Convolve the neural model by the BOLD HRF model

yBOLD = conv(yDN,BOLDHRF);
yBOLD = yBOLD(1:modelLength); % NEED TO DETERMINE WHY THIS IS NECESSARY

% Plot the model stages

figure
hold on;
plot(neuralIRF);
plot(decayingExponential);
plot(BOLDHRF);
plot(BOLDHRF,'.');
hold off;

figure;
hold on;
plot(yStimulus);
plot(yLinear);
plot(yCTS);
plot(yDN);
plot(yBOLD);
hold off;
