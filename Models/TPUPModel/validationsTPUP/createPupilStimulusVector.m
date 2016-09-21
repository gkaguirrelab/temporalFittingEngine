function [stimulus] = createPupilStimulusVector()

% parameters of the stimulus simulation
modelDuration=13; % total duration of the modeled period, in seconds
modelResolution=1000; % temporal sampling frequency of model, in Hz
plateauDuration=2; % duration of the modeled step of neural activity, in seconds
rampDuration=0.5; % duration of the half-cosine ramp on and off

% Create a stimulus model, which is a simple step function of input,
% with a half cosine ramp on and off,
t = linspace(0,modelDuration,modelDuration*modelResolution); % the temporal domain of the model
halfCosineOff=(cos(linspace(0,pi,rampDuration*modelResolution))+1)/2;
halfCosineOn=fliplr(halfCosineOff);

rampLength=length(halfCosineOn);
plateauLength=plateauDuration*modelResolution;

% Build the stimulus vector
stimulus = t*0; % the stimulus vector
stimulus(1:rampLength)=halfCosineOn;
stimulus(1+rampLength:rampLength+plateauLength)=1;
stimulus(1+rampLength+plateauLength:rampLength*2+plateauLength)=halfCosineOff;

% return column vectors
stimulus=stimulus';

end

