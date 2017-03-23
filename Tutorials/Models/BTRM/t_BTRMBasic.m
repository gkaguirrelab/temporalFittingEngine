function [ paramsFit ] = t_BTRMBasic(varargin)
% function [ paramsFit ] = t_BTRMBasic(varargin)
%
% Demonstrate function for the Block Temporal Response Model model.
%
% The model inspired by:
%   Zaidi, Q., Ennis, R., Cao, D., & Lee, B. (2012). Neural locus of color afterimages. Current Biology, 22(3), 220-224.
%   Zhou, J., Benson, N. C., Kay, K., & Winawer, J. (2017). Systematic changes in temporal summation across human visual cortex. bioRxiv, 108639.
%
% Optional key/value pairs
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});

%% Demo 1 -- figure 3A of Zaidi
fprintf('Demo 1 -- Figure 3A of Zaidi - Half-sine stimulus subjected to a inhibitory feedback with 8 second time constant.\n');
fprintf('  In blue is the stimulus, in grey the simulated response, in red the model fit.\n\n');

% Housekeeping and setup
temporalFit = tfeBTRM('verbosity','none'); % Construct the model object

% Create a half-sinusoid stimulus
deltaT = 100; % in msecs
totalTime = 34000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
stimulusStruct.values = sin(2*pi*stimulusStruct.timebase / totalTime);
stimulusStruct.values(stimulusStruct.values<0)=0;
defaultParamsInfo.nInstances = 1;

% Model params will remove the dCTS, minimize the CTS components, and
% select other parameters similar to Zaidi
params=temporalFit.defaultParams('defaultParamsInfo',defaultParamsInfo,'use_dCTS',false);
params.paramMainMatrix(1)=50; % spikes per second
params.paramMainMatrix(2)=1; % minimal gamma IRF effect (this is an impulse)
params.paramMainMatrix(3)=8; % 8 second time constant of leaky negative integrator
params.paramMainMatrix(4)=0.25; % the inhibitory component is 25% of the positive effect
params.paramMainMatrix(5)=1; % no compressive non-linearity effect

params.noiseSd=5; % stdev of noise
params.noiseInverseFrequencyPower=0; % white noise

if p.Results.generatePlots
    % Create a figure window
    figure;
    % plot the stimulus profile and a refline
    plot(stimulusStruct.timebase/1000,stimulusStruct.values(1,:)*params.paramMainMatrix(1),'-b','DisplayName','stimulus');
    hold on
    refline(0,0);
end

% create the simulated response and plot it
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,[],'AddNoise',true);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 0.5 0.5],'DisplayName','spikes per second');
end

% Construct a packet for fitting
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = []; thePacket.metaData = [];

% Fit the simulated data and plot the result
[paramsFit.demo1,fVal,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'diffMinChange',0.01,...
    'defaultParamsInfo', defaultParamsInfo);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[1 0.25 0.25],'DisplayName','fit');
end


%% Demo 2 -- figure 4A of Zhou
fprintf('Demo 2 -- Figure 4A of Zhou - stimuli (black) with log spaced (x2) durations.\n');
fprintf('  In green is the BOLD response to a linear model, in blue with a compressive non-linearity.\n\n');

% Housekeeping and setup
stimulusStruct=[];
thePacket=[];
modelResponseStruct=[];
params=[];

temporalFit = tfeBTRM('verbosity','none'); % Construct the model object

% Create a set of step stimuli with varying durations
stimDurations=[17,33,67,134,267,533];
deltaT = 1; % in msecs
totalTime = 16000*(length(stimDurations)+1); % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
for ii=1:length(stimDurations)
    stimulusStruct.values(ii,:)=zeros(1,length(stimulusStruct.timebase));
    stimulusStruct.values(ii,16000*ii+1:16000*ii+stimDurations(ii))=1;
end
defaultParamsInfo.nInstances = length(stimDurations);

% Model params will remove the dCTS, use the CTS components, and
% minimize Zaidi
params=temporalFit.defaultParams('defaultParamsInfo',defaultParamsInfo,'use_dCTS',false);
params.paramMainMatrix(:,1)=1; % neural units of response (e.g., % change)
params.paramMainMatrix(:,2)=75; % gamma IRF effect measured for V1
params.paramMainMatrix(:,3)=8; % 8 second time constant of leaky negative integrator
params.paramMainMatrix(:,4)=0; % no inhibitory delayed component
params.paramMainMatrix(:,5)=0.27; % compressive non-linearity observed in V1

params.noiseSd=0.0001; % stdev of noise
params.noiseInverseFrequencyPower=1; % pink noise

% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
kernelStruct.timebase=stimulusStruct.timebase;
% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;
% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct.values=kernelStruct.values-kernelStruct.values(1);
kernelStruct=normalizeKernelArea(kernelStruct);

if p.Results.generatePlots
    % Create a figure window
    figure;
    hold on
    % plot the stimulus profiles and a refline
    for ii=1:length(stimDurations)
        plot(stimulusStruct.timebase/1000,stimulusStruct.values(ii,:)*params.paramMainMatrix(1)/5,'-k','DisplayName','stimulus');
    end
    refline(0,0);
end

% create the simulated response and plot it, with a compressive
% non-linearity
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,kernelStruct,'AddNoise',true);
if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 0.5 1],'DisplayName','with compressive nonlinearity');
end

% Do it again without the non-linearity
params.paramMainMatrix(:,5)=1;
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,kernelStruct,'AddNoise',true);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 1 0.5],'DisplayName','w/o compressive nonlinearity');
end

%% Demo 3 -- figure 6A of Zhou
fprintf('Demo 3 -- Figure 6A of Zhou - Step function stimulus subjected to the dCTS model with params from V1.\n');
fprintf('  In blue is the stimulus, in grey the simulated response, in red the model fit.\n\n');

% Housekeeping and setup
stimulusStruct=[];
thePacket=[];
modelResponseStruct=[];
params=[];

% Create a 500 msec step stimulus
deltaT = 1; % in msecs
totalTime = 1250; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
stimulusStruct.values = zeros(1,length(stimulusStruct.timebase));
stimulusStruct.values(1,250:750)=1;
defaultParamsInfo.nInstances = 1;

% We use the dCTS, and minimize the Zaidi component
params=temporalFit.defaultParams('defaultParamsInfo',defaultParamsInfo,'use_dCTS',true);
params.paramMainMatrix(1)=15; % broadband power
params.paramMainMatrix(2)=90; % gamma IRF time constant in msecs
params.paramMainMatrix(3)=8; % 8 second time constant of leaky negative integrator
params.paramMainMatrix(4)=0; % the inhibitory component is 0% of the positive effect
params.paramMainMatrix(5)=1.8; % compression in the dCTS
params.paramMainMatrix(6)=0.1; % adaptive time constant (in seconds)
params.paramMainMatrix(7)=0.1; % sigma saturation constant

params.noiseSd=1; % stdev of noise
params.noiseInverseFrequencyPower=0; % white noise

if p.Results.generatePlots
    % Create a figure window
    figure;
    % plot the stimulus profile and a refline
    plot(stimulusStruct.timebase/1000,stimulusStruct.values(1,:)*params.paramMainMatrix(1),'-b','DisplayName','stimulus');
    hold on
    refline(0,0);
end

% create the simulated response and plot it
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,[],'AddNoise',true,'use_dCTS',true);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 0.5 0.5],'DisplayName','broad band power');
end

% Construct a packet for fitting
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = []; thePacket.metaData = [];

% Fit the simulated data and plot the result
[paramsFit.demo3,fVal,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo,...
    'use_dCTS',true);
if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[1 0.25 0.25],'DisplayName','fit');
end


%% Demo 4 -- random series of pulses
fprintf('Demo 4 -- 3 second stimulus pulses (light blue) with in random binary sequence.\n');
fprintf('  In blue/gray is the prediction of the CTS model with long-term inibition (neural and BOLD), in red is the model fit.\n\n');

% Housekeeping and setup
stimulusStruct=[];
thePacket=[];
modelResponseStruct=[];
params=[];

% Define an event that is a single, windowed pulse
deltaT = 100;
eventDuration=3000; % msecs
rampDuration=500; % msecs
% the square wave step
eventStruct.values=zeros(1,eventDuration/deltaT);
eventStruct.values(1:round(eventDuration/deltaT))=1;
% half cosine ramp on
eventStruct.values(1:round(rampDuration/deltaT))= ...
    (fliplr(cos(linspace(0,2*pi,round(rampDuration/deltaT))/2))+1)/2;
% half cosine ramp off
eventStruct.values(1+round(eventDuration/deltaT)-round(rampDuration/deltaT): ...
    round(eventDuration/deltaT))= ...
    (cos(linspace(0,2*pi,round(rampDuration/deltaT))/2)+1)/2;

% Build a binary array of these events
totalTime = 300000; % in msecs. This is a 5 minute experiment
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);
stimulusStruct.values = zeros(1,nTimeSamples);
nEvents=totalTime/eventDuration;
for ii=1:nEvents
    if rand()>0.5
        stimulusStruct.values(1,eventDuration*(ii-1)/deltaT+1:(eventDuration*(ii)/deltaT))=eventStruct.values;
    end
end
defaultParamsInfo.nInstances = 1;

% Model params will remove the dCTS, but use CTS and Zaidi
params=temporalFit.defaultParams('defaultParamsInfo',defaultParamsInfo,'use_dCTS',false);
params.paramMainMatrix(:,1)=1; % neural units of response (e.g., % change)
params.paramMainMatrix(:,2)=75; % gamma IRF effect measured for V1
params.paramMainMatrix(:,3)=12; % 12 second time constant of leaky negative integrator
params.paramMainMatrix(:,4)=0.25; % inhibitory delayed component equal to 25% of primary response
params.paramMainMatrix(:,5)=0.27; % compressive non-linearity observed in V1

params.noiseSd=0.05; % stdev of noise
params.noiseInverseFrequencyPower=1; % pink noise

% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
kernelStruct.timebase=stimulusStruct.timebase;
% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;
% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct.values=kernelStruct.values-kernelStruct.values(1);
kernelStruct=normalizeKernelArea(kernelStruct);

if p.Results.generatePlots
    % Create a figure window
    figure;
    % plot the stimulus profiles and a refline
    plot(stimulusStruct.timebase/1000,stimulusStruct.values(1,:)*params.paramMainMatrix(1),'-','Color',[0.75 0.75 1],'DisplayName','stimulus');
    hold on
    refline(0,0);
end

% create the simulated without the HRF kernel or noise response and plot it
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,[],'AddNoise',false);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 0.5 1],'DisplayName','simulated response');
end

% create the simulated response and plot it
modelResponseStruct = temporalFit.computeResponse(params,stimulusStruct,kernelStruct,'AddNoise',true);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[0.5 0.5 0.5],'DisplayName','simulated response');
end

% Construct a packet for fitting
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];

% Fit the simulated data and plot the result
[paramsFit.demo4,fVal,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo,...
    'use_dCTS',false);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'NewWindow',false,'Color',[1 0.25 0.25],'DisplayName','fit');
end

% Report the modeled and observed parameters
fprintf('Simulated parameters:\n');
temporalFit.paramPrint(params);
fprintf('\n');
fprintf('Measured parameters:\n');
temporalFit.paramPrint(paramsFit.demo4);
fprintf('\n');


end % function