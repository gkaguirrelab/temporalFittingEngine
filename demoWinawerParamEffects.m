%% Demo the effects of parameter choices in Winawer model
%

% parameters of the stimulus simulation
modelDuration=40; % total duration of the modeled period, in seconds
modelResolution=1; % temporal sampling frequency of model, in Hz
stimulusDuration=24; % duration of the modeled step of neural activity, in seconds

% Create a stimulus model, which is a simple step function of input
t = linspace(1,modelDuration,modelDuration*modelResolution); % the temporal domain of the model
yStimulus = t*0; % the stimulus vector
yStimulus(2:(stimulusDuration*modelResolution)+1) = 1; % implement the step function of input

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

% Show the effect of varying sigma, n, tau2

figure;

% vary sigma
holdParam=param.sigma;
subplot(1,3,1);
hold on;
for i=0.1:.1:.5
    param.sigma=i;
    yBOLD=fitWinawerModelToBOLD(t,yStimulus,false,param);
    plot(yBOLD);
end
title('Vary sigma, 0.1:0.1:0.5');
hold off;
param.sigma=holdParam;

% vary n
holdParam=param.n;
subplot(1,3,2);
hold on;
for i=0.5:.1:1
    param.n=i;
    yBOLD=fitWinawerModelToBOLD(t,yStimulus,false,param);
    plot(yBOLD);
end
title('Vary n, 0.5:0.1:1');
hold off;
param.n=holdParam;

% vary tau2
holdParam=param.tau2;
subplot(1,3,3);
hold on;
for i=-2:1:0
    param.tau2=10^i;
    yBOLD=fitWinawerModelToBOLD(t,yStimulus,false,param);
    plot(yBOLD);
end
title('Vary tau2, 0.01, 0.1, 1');
hold off;
param.tau2=holdParam;