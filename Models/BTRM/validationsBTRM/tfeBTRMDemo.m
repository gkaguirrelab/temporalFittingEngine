% tfeBlockTemporalResponseMRIModelDemo
%
% Demonstrate function for the BTRM Model.
%

%% Clear and close
clearvars;
%close all;

% Reset the random number generator for consistent results
%rng default

%% Construct the model object
temporalFit = tfeBTRM('verbosity','none');

%% Get the default forward model parameters
params0 = temporalFit.defaultParams();
fprintf('Default model parameters:\n');
temporalFit.paramPrint(params0);
fprintf('\n');

%% Temporal domain of the stimulus 
deltaT = 100; % in msecs
totalTime = 50000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify a stimulus. 
% We create here a step function of neural activity, with half-cosine ramps
%  on and off
eventOnset=0;
eventDuration=3000; % msecs
rampDuration=500; % msecs

% the square wave step
eventStruct.values=zeros(1,eventDuration/deltaT);
eventStruct.values(round(eventOnset/deltaT)+1: ...
                      round(eventOnset/deltaT)+round(eventDuration/deltaT))=1;
% half cosine ramp on
eventStruct.values(round(eventOnset/deltaT)+1: ...
                      round(eventOnset/deltaT)+round(rampDuration/deltaT))= ...
                      (fliplr(cos(linspace(0,2*pi,round(rampDuration/deltaT))/2))+1)/2;
% half cosine ramp off
eventStruct.values(round(eventOnset/deltaT)+1+round(eventDuration/deltaT)-round(rampDuration/deltaT): ...
                      round(eventOnset/deltaT)+round(eventDuration/deltaT))= ...
                      (cos(linspace(0,2*pi,round(rampDuration/deltaT))/2)+1)/2;

% This stimulus is just the positive portion of a sinusoid with a cycle
% time equal to totalTime.
% stimulusStruct.values = sin(2*pi*stimulusStruct.timebase / totalTime);
% stimulusStruct.values(stimulusStruct.values<0)=0;


%% Temporal domain of the stimulus 
deltaT = 100; % in msecs
totalTime = 600000; % in msecs. This is a 10 minute experiment
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);
stimulusStruct.values = zeros(1,nTimeSamples);

nEvents=totalTime/eventDuration;
for ii=1:nEvents
    if rand()>0.5
    stimulusStruct.values(1,eventDuration*(ii-1)/deltaT+1:(eventDuration*(ii)/deltaT))=eventStruct.values;
    end
end

%% Define a kernelStruct. In this case, a double gamma HRF
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

%% Create and plot modeled responses

% Set the noise level and report the params
params0.noiseSd = 0.2;
fprintf('Simulated model parameters:\n');
temporalFit.paramPrint(params0);
fprintf('\n');

% Create a figure window
figure;

% First create and plot the response without noise and without convolution
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,[],'AddNoise',false);
temporalFit.plot(modelResponseStruct,'NewWindow',false,'DisplayName','neural response','Color',[.5 .5 1]);
hold on;

% Add the stimulus profile to the plot
plot(stimulusStruct.timebase/1000,stimulusStruct.values(1,:),'-k','DisplayName','stimulus');

% Now plot the response with convolution and noise, as well as the kernel
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',true);
temporalFit.plot(modelResponseStruct,'NewWindow',false,'DisplayName','noisy BOLD response');
plot(kernelStruct.timebase/1000,kernelStruct.values/max(kernelStruct.values),'-b','DisplayName','kernel');

%% Construct a packet and model params
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;

%% Test the fitter
[paramsFit,fVal,modelResponseStruct] = ...
            temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo);
        
%% Report the output
fprintf('Model parameter from fits:\n');
temporalFit.paramPrint(paramsFit);
fprintf('\n');

temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit');
legend('show');legend('boxoff');
hold off;

