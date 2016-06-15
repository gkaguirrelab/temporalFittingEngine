


% parameters of the stimulus simulation

modelDuration=336; % total duration of the modeled period, in seconds
modelResolution=10; % temporal sampling frequency of model, in Hz
modelLength=modelDuration*modelResolution; % duh
stimulusDuration=12; % duration of the modeled step of neural activity, in seconds
stimulusOnsets=[1,36,112,124,172,213,306]; % the time of onset of stimulus blocks, in seconds


% Create a stimulus model, which is a simple step function of input
t = linspace(1,modelDuration,modelDuration*modelResolution); % the temporal domain of the model
yStimulus = t*0; % the stimulus vector
for i=1:length(stimulusOnsets)
    thisStimulusOnset=round(stimulusOnsets(i)*modelResolution);    
    yStimulus(thisStimulusOnset:thisStimulusOnset+(stimulusDuration*modelResolution)+1) = 1; % implement the step function of input
end

% parameters of the double-gamma hemodynamic filter (HRF)
param.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in seconds)
param.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in seconds)
param.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

%% Convolve the neural model by the BOLD HRF

% Create a double-gamma model of the hemodynamic response function (HRF)
BOLDHRF = gampdf(t, param.gamma1, 1) - ...
    gampdf(t, param.gamma2, 1)/param.gammaScale;

% scale to unit sum to preserve amplitude of y following convolution
BOLDHRF = BOLDHRF/sum(BOLDHRF);

% Perform the convolution
yBOLD = conv(yStimulus,BOLDHRF);

% Truncate the convolved vector to the input length. Not sure why I have to
% do this. Probably mis-using the conv function.
yBOLD = yBOLD(1:modelLength);



%% Plot the model
figure;
hold on;
r1 = plot(yStimulus);
r2 = plot(yBOLD);
title('Stimulus and HRF convolved vector');
legend([r1 r2], 'stimulus', 'HRF_convolved');
legend boxoff
hold off;