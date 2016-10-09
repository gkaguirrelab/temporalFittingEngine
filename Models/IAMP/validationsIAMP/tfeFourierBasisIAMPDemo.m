% tfeInstanceAmplitudeDemo
%
% Demonstrate the creation and application of a Fourier basis set analysis.
%

%% Clear and close
clear; close all;

%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

%% Temporal domain of the stimulus model
totalTime = 330000; % in msecs. This is a 5:30 duration experiment
stimulusStruct.timebase = linspace(0,totalTime-1,totalTime);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus 
% We will create a set of impulses of various amplitudes in a stimulus
% matrix. They will be jittered in timing ±5 seconds, uniformly distributed.
eventTimes=linspace(6000,306000,21);
eventJitter=round((rand(1,21)-.5)*15000);
eventTimes=eventTimes+eventJitter;
eventTimes(eventTimes<0)=0;
eventTimes(eventTimes>(totalTime-1000))=0;
nInstances=length(eventTimes);
eventDuration=1; % pulse duration in msecs
defaultParamsInfo.nInstances = nInstances;
stimulusStruct.values=stimulusStruct.timebase*0;
for ii=1:nInstances
    stimulusStruct.values(1,eventTimes(ii):eventTimes(ii)+eventDuration)=1;
end

%% Define a double gamma HRF kernel struct
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
kernelStruct.timebase =linspace(0,15999,16000);
kernelStruct.values = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct=prepareHRFKernel(kernelStruct);

%% Get the default forward model parameters
params0 = temporalFit.defaultParams('defaultParamsInfo', defaultParamsInfo);
%params0.noiseSd = 0.000002;

%% Create a simulated response
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',false);
newTimebase=linspace(0,330000-1000,330000/1000);
modelResponseStruct = temporalFit.resampleTimebase(modelResponseStruct,newTimebase);
modelResponseStruct.values=modelResponseStruct.values-mean(modelResponseStruct.values);

%% Plot the simulated response
figure;
temporalFit.plot(modelResponseStruct,'NewWindow',false,'DisplayName','neural response','Marker','o');
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is where we build the Fourier model and apply it to a modeled
% response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Create a fourier StimStructure
[ stimulusStruct, fourierSetStructure ] = ...
    makeFourierStimStruct( stimulusStruct.timebase, eventTimes, 16000, 16 );

%% Construct a packet and model params
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Get the default forward model parameters
params0 = temporalFit.defaultParams('defaultParamsInfo', defaultParamsInfo);

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We treat each Fourier component as an instance
defaultParamsInfo.nInstances = 16;

%% Test the fitter
[paramsFit,fVal,modelResponseStruct] = ...
            temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix, ...
            'searchMethod','linearRegression');

%% Make some plots
% Plot of the temporal fit results
temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit');
legend('show');legend('boxoff');
hold off

% Obtain the high-resolution estimate of the kernel
estimatedKernelStruct.timebase = kernelStruct.timebase;
estimatedKernelStruct.values = ...
    (fourierSetStructure.values'*paramsFit.paramMainMatrix)';
estimatedKernelStruct = prepareHRFKernel(estimatedKernelStruct);

% Plot the true kernel and the estimated kernel.
% These are not expected to be identical in this simulation as the
% high-resolution information is available through a finite number of
% jittered events. Try turning off event jitter to see how the estimate is
% made worse.
figure
plot(estimatedKernelStruct.timebase,estimatedKernelStruct.values,'Color',[0 1 0],'DisplayName','estimate')
hold on
plot(kernelStruct.timebase,kernelStruct.values,'Color',[1 0 0],'DisplayName','true')
legend('show');legend('boxoff');
hold off
