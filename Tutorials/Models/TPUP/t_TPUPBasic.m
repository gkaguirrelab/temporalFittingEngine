function validationData = t_TPUPBasic(varargin)
% validationData = t_TPUPBasic(varargin)
%
% Demonstrate the Two Component PUPil response model
%
% Optional key/value pairs
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});


%% Construct the model object
temporalFit = tfeTPUP('verbosity','none');

%% Get the default forward model parameters
params0 = temporalFit.defaultParams;
fprintf('Default model parameters:\n');
temporalFit.paramPrint(params0);
fprintf('\n');

%% Temporal domain of the stimulus
deltaT = 10; % in msecs
totalTime = 14000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus struct.
% We create here a step function of neural activity, with half-cosine ramps
%  on and off
stepOnset=1000; % msecs
stepDuration=3000; % msecs
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

%% We will not make use of a kernel in this model
kernelStruct=[];

%% Create a modeled pupil response, with added noise
params0.noiseSd = 0.05;
fprintf('Simulated model parameters:\n');
temporalFit.paramPrint(params0);
fprintf('\n');

modelResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',true);

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct);
    hold on
    plot(stimulusStruct.timebase/1000,stimulusStruct.values(1,:),'-k');
end

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

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false);
    hold off
end

%% Set returned validationData structure
if (nargout > 0)
    validationData.params1 = paramsFit;
    validationData.modelResponseStruct = modelResponseStruct;
    validationData.thePacket = thePacket;
end

end
