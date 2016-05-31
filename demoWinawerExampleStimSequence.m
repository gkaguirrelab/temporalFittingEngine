%% Demo a hypothetical stimulus sequence with variation in temporal parameters.
%
%  05-31-2016 -- gka wrote it

% parameters of the stimulus simulation
modelResolution=1; % temporal sampling frequency of model, in Hz
stimulusDuration=12; % duration of the modeled step of neural activity, in seconds

% parameters of the neural filters
param.tau1 = 0.005;   % time constant of the neural IRF (in seconds)
param.epsilon = 0.35; % compressive non-linearity parameter
param.tau2 = 0.1;     % time constant of the low-pass (exponential decay) component of delayed normalization
param.sigma = 0.13;   % constant scaling factor of the divisive normalization
param.n = .5;          % the order of the delayed normalization

% parameters of the double-gamma hemodynamic filter (HRF)
param.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in seconds)
param.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in seconds)
param.gammaScale = 3; % scaling factor between the positive and negative gamma componenets

% a vector of amplitudes of neural response for each stimulus block

simulationGain=[1,1,0,1,1,0,1,1,0,1,1];
simulationTau2=[0.1,0.1,0,0.1,1,0,1,0.1,0,1,1];
    
modelDuration=stimulusDuration*length(simulationAmplitudes); % total duration of the modeled period, in seconds
t = linspace(1,modelDuration,modelDuration*modelResolution); % the temporal domain of the model
fullModel=t*0;

figure
hold on;
for i=1:length(simulationAmplitudes);
yStimulus = t*0; % the stimulus vector
yStimulus((i-1)*(stimulusDuration*modelResolution)+1:i*(stimulusDuration*modelResolution)) = 1; % implement the step function of input
param.tau2=simulationTau2(i);
yBOLD=fitWinawerModelToBOLD(t,yStimulus,false,param);
yBOLD=yBOLD*simulationGain(i);
plot(yBOLD);
fullModel=fullModel+yBOLD;
end
hold off;

% Show the resulting simulated signal

figure;
plot(fullModel);
