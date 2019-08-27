function validationData = t_tfeRegressWithNan(varargin)
% validationData = t_tfeRegressWithNan(varargin)
%
% Demonstrate the performance of the linear regression operation when the
% response vector contains nan values
%


%% Parse vargin for options passed here
p = inputParser; p.PartialMatching = false;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});


%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

%% Temporal domain of the stimulus
deltaT = 100; % in msecs
totalTime = 330000; % in msecs. This is a 5:30 duration experiment
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);
eventDuration=10000; % pulse duration in msecs
nInstances = 1;

%% Specify the stimulus struct.
% This is a simple square wave of stimulation every 20 seconds
eventTimes=deltaT:eventDuration*2:floor(totalTime/(eventDuration*2))*(eventDuration*2)+deltaT;
stimulusStruct.values(1,:)=zeros(1,nTimeSamples);
for ii=1:length(eventTimes)
    stimulusStruct.values(1,eventTimes(ii)/deltaT:(eventTimes(ii)/deltaT)+eventDuration/deltaT-1)=1;
end
defaultParamsInfo.nInstances = nInstances;


%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

kernelStruct.timebase=linspace(0,15999,16000);

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.  It's defined on 1 msec sampling
% timebase.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% Normalize the kernel to have unit amplitude
[ kernelStruct ] = normalizeKernelArea( kernelStruct );

% When the IAMP model computes the response, it resamples the kernal to the
% same timebase as the response, if they differ.  That's slow.  So we can
% speed things up by doing it once here.
if (deltaT ~= 1)
    nSamples = ceil((kernelStruct.timebase(end)-kernelStruct.timebase(1))/deltaT);
    newKernelTimebase = kernelStruct.timebase(1):deltaT:(kernelStruct.timebase(1)+nSamples*deltaT);
    kernelStruct = temporalFit.resampleTimebase(kernelStruct,newKernelTimebase);
end

%% Get the default forward model parameters
params0 = temporalFit.defaultParams('defaultParamsInfo', defaultParamsInfo);

%% Create a modeled response with noise
% Set the noise level and report the params
params0.noiseSd = 0.02;
modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',true);


%% Construct a packet and model params
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instamce.
defaultParamsInfo.nInstances = nInstances;

%% Obtain regression values with the intact response vector
[paramsFit,~,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, ...
    'searchMethod','linearRegression');

%% Add nans to 10% of the response values at random locations
nanIdx = randsample(length(thePacket.response.values),floor(length(thePacket.response.values)*0.1));
thePacket.response.values(nanIdx) = nan;

%% Obtain regression values with the nan-ed response vector
[paramsFitNan,~,modelResponseStructNan] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, ...
    'searchMethod','linearRegression');


% Plot of the temporal fit results
if p.Results.generatePlots
    temporalFit.plot(thePacket.response,'Color',[0 0 1],'NewWindow',true,'DisplayName','response');
    temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit without nans');
    temporalFit.plot(modelResponseStructNan,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit with nans');
    legend('show');legend('boxoff');
    hold off
end

%% Set validation data for return
if (nargout > 0)
    validationData.responseStruct = responseStruct;
    validationData.kernelStruct = kernelStruct;
    validationData.convResponseStruct = convResponseStruct;
    validationData.resampledKernelStruct = resampledKernelStruct;
end

end

