% tfeBlockTemporalResponseMRIModelDemo
%
% Demonstrate function for the BTRM Model.
%

%% Clear and close
clear; close all;

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

%% Specify the stimulus values. 
% We create here a step function of neural activity, with half-cosine ramps
%  on and off
stepOnset=1000; % msecs
stepDuration=12000; % msecs
rampDuration=3000; % msecs

% the square wave step
stimulusStruct.values=zeros(1,nTimeSamples);
stimulusStruct.values(round(stepOnset/deltaT): ...
                      round(stepOnset/deltaT)+round(stepDuration/deltaT)-1)=1;
% half cosine ramp on
stimulusStruct.values(round(stepOnset/deltaT): ...
                      round(stepOnset/deltaT)+round(rampDuration/deltaT)-1)= ...
                      fliplr(cos(linspace(0,pi,round(rampDuration/deltaT))/2));
% half cosine ramp off
stimulusStruct.values(round(stepOnset/deltaT)+round(stepDuration/deltaT)-round(rampDuration/deltaT): ...
                      round(stepOnset/deltaT)+round(stepDuration/deltaT)-1)= ...
                      cos(linspace(0,pi,round(rampDuration/deltaT))/2);

%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

kernelStruct.timebase=stimulusStruct.timebase;

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;

% scale to unit sum to preserve amplitude of signal following convolution
hrf=hrf/sum(hrf);
kernelStruct.values=hrf;

%% Create and plot modeled responses

% Set the noise level and report the params
params0.noiseSd = 0.05;
fprintf('Simulated model parameters:\n');
temporalFit.paramPrint(params0);
fprintf('\n');

% Create a figure window
figure;

% First create and plot the response without noise and without convolution
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,[],'AddNoise',false);
temporalFit.plot(modelResponseStruct,'NewWindow',false,'DisplayName','neural response');
hold on

% Add the stimulus profile to the plot
plot(stimulusStruct.timebase,stimulusStruct.values(1,:),'-k','DisplayName','stimulus');

% Now plot the response with convolution and noise, as well as the kernel
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',true);
temporalFit.plot(modelResponseStruct,'NewWindow',false,'DisplayName','noisy BOLD response');
plot(kernelStruct.timebase,kernelStruct.values/max(kernelStruct.values),'-b','DisplayName','kernel');

%% Construct a packet and model params
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instamce.
defaultParamsInfo.nInstances = 1;

%% Test the fitter
[paramsFit,fVal,modelResponseStruct] = ...
            temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix);
        
%% Report the output
fprintf('Model parameter from fits:\n');
temporalFit.paramPrint(paramsFit);
fprintf('\n');

temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit');
legend('show');legend('boxoff');
hold off

