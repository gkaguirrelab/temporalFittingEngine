% tfeBlockTemporalResponseMRIModelDemo
%
% Demonstrate function for the BTRM Model.
%

%% Clear and close
clear; close all;

%% Construct the model object
temporalFit = tfeBTRM('verbosity','high');

%% Get the default forward model parameters
params0 = temporalFit.defaultParams;
fprintf('Default model parameters:\n');
temporalFit.paramPrint(params0);

%% Temporal domain of the stimulus 
deltaT = 100; % in msecs
totalTime = 26000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus struct. 
% We create here a step function of neural activity, with half-cosine ramps
%  on and off
stepOnset=1000; % msecs
stepDuration=12000; % msecs
rampDuration=500; % msecs

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

hrfTimebaseSecs=stimulusStruct.timebase/1000; % in seconds
hrf = gampdf(hrfTimebaseSecs+1, hrfParams.gamma1, 1) - ...
    gampdf(hrfTimebaseSecs+1, hrfParams.gamma2, 1)/hrfParams.gammaScale;

% scale to unit sum to preserve amplitude of signal following convolution
hrf=hrf/sum(hrf);
kernelStruct.values=hrf;
kernelStruct.timebase=stimulusStruct.timebase;

%% Create a modeled fMRI response, with added noise
params0.noiseSd = 0.02;
fprintf('Simulated model parameters:\n');
temporalFit.paramPrint(params0);
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',false);
temporalFit.plot(modelResponseStruct);
hold on
plot(stimulusStruct.timebase,stimulusStruct.values(1,:),'-k');
plot(kernelStruct.timebase,kernelStruct.values/max(kernelStruct.values),'-b');

% 
% 
% %% Construct a packet
% thePacket.stimulus = stimulusStruct;
% thePacket.response = modelResponseStruct;
% thePacket.kernel = [];
% thePacket.metaData = [];
% 
% %% Test the fitter
% [paramsFit,fVal,fitResponseStruct] = temporalFit.fitResponse(thePacket);
% fprintf('Model parameter from fits:\n');
% temporalFit.paramPrint(paramsFit);
% temporalFit.plot(fitResponseStruct,'Color',[0 1 0],'NewWindow',false);

