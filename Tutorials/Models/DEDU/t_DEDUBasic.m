function validationData = t_DEDUBasic(varargin)
% validationData = t_DEDUBasic(varargin)
%
% Demonstrate the DElay and DUration model
%
% We will model a single event.
%
% Optional key/value pairs
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});

%% Construct the model object
tfeHandle = tfeDEDU('verbosity','none');

%% Temporal domain of the stimulus
deltaT = 1; % in msecs
totalTime = 20000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

% The stimulus profile is just a delta function. The DEDU model discards
% the stimulus profile and simply identifies the time of onset of the
% stimulus and uses this.
nInstances=1;
defaultParamsInfo.nInstances=nInstances;
stimulusStruct.values(1,:)=zeros(1,nTimeSamples);
stimulusStruct.values(1,1)=1;


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

% Get the default forward model parameters
params0 = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
params0.noiseSd = 0.01;

% start the packet assembly
thePacket.stimulus = stimulusStruct;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];

% Randomize the order of the stimTypes
thePacket.stimulus.metaData.stimTypes=[1];
thePacket.stimulus.metaData.stimLabels=['demo'];

% Create some params to define the simulated data for this packet
paramsLocal=params0;
paramsLocal.paramMainMatrix(1,1)=1;
paramsLocal.paramMainMatrix(1,2)=0.2; % delay of onset of response in seconds
paramsLocal.paramMainMatrix(1,3)=3.5; % duration of response in seconds

%% Report the modeled params
fprintf('Simulated model parameters:\n');
tfeHandle.paramPrint(paramsLocal);
fprintf('\n');

% Generate the simulated response
simulatedResponseStruct = tfeHandle.computeResponse(paramsLocal,thePacket.stimulus,thePacket.kernel,'AddNoise',true);

% Add the simulated response to this packet
thePacket.response=simulatedResponseStruct;

if p.Results.generatePlots
    tfeHandle.plot(simulatedResponseStruct,'DisplayName','Simulated');
end

%% Test the fitter
% The DEDU model requires a diff min change, as the duration parameter is
% discretized to have effects in untits of msecs. Failure to do so will
% lead fmincon to believe that changes in the duration parameter have no
% effect upon the model fit.
[paramsFit,fVal,modelResponseStruct] = ...
    tfeHandle.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, ...
    'DiffMinChange',0.001);

%% Report the output
fprintf('Model parameter from fits:\n');
tfeHandle.paramPrint(paramsFit);
fprintf('\n');

if p.Results.generatePlots
    tfeHandle.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit');
    legend('show');legend('boxoff');
end

%% Set returned validationData structure
if (nargout > 0)
    validationData.params1 = paramsFit;
    validationData.modelResponseStruct = modelResponseStruct;
    validationData.thePacket = thePacket;
end

end % function
