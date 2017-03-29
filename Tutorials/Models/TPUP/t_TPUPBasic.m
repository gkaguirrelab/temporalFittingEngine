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

%% Supply some model parameters that make a nice cone-driven response:
params0.paramMainMatrix=[151.8345, 183.3596, 3.0000, -5.3529, -31.1735, -40.5347];

%% Temporal domain of the stimulus
deltaT = 1; % in msecs
totalTime = 13000; % in msecs
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
    fliplr((cos(linspace(0,pi*2,round(rampDuration/deltaT))/2)+1)/2);
% half cosine ramp off
stimulusStruct.values(round(stepOnset/deltaT)+round(stepDuration/deltaT)-round(rampDuration/deltaT): ...
    round(stepOnset/deltaT)+round(stepDuration/deltaT)-1)= ...
    (cos(linspace(0,pi*2,round(rampDuration/deltaT))/2)+1)/2;

%% We will not make use of a kernel in this model
kernelStruct=[];

%% Create a modeled pupil response, with added noise
params0.noiseSd = .5;
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


% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instamce.
defaultParamsInfo.nInstances = 1;

%% Test the fitter
[paramsFit,fVal,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo);

%% Report the output
fprintf('Model parameter from fits:\n');
temporalFit.paramPrint(paramsFit);
fprintf('\n');

if p.Results.generatePlots
    temporalFit.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false);
    hold off
end

%% Create a figure that illustrates the model components
if p.Results.generatePlots
    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    % Plot the stimulus and the stimulusSlewOn
    h=subplot(2,2,1);
    temporalFit.plot(stimulusStruct,'Color',[.5 .5 .5],'NewWindow',false,'DisplayName', 'stimulus');
    pbaspect([1 2 1])
    xlim([0 13]);
    set(h,'Xtick',0:2:13)
    hold on
    tmpParams=params0;
    tmpParams.paramMainMatrix([5,6])=0;
    tmpParams.paramMainMatrix([4,2])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[.75 .75 .75],'NewWindow',false,'DisplayName', 'slewRateOn');
    legend('show')
    hold off
    % Plot the gamma and exponential kernels
    h=subplot(2,2,2);
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,6])=0;
    tmpParams.paramMainMatrix([3,5])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,makeImpulseStimStruct(stimulusStruct),[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 0 0],'NewWindow',false,'DisplayName', 'gamma');
    pbaspect([1 2 1])
    xlim([0 13]);
    set(h,'Xtick',0:2:13)
    hold on
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,5])=0;
    tmpParams.paramMainMatrix([2,6])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,makeImpulseStimStruct(stimulusStruct),[]);
    temporalFit.plot(modelResponseStruct,'Color',[.5 0 0],'NewWindow',false,'DisplayName', 'exponential');
    legend('show')
    hold off
    % Plot the transient, sustained, and persistent components with unit
    % area
    h=subplot(2,2,3);
    tmpParams=params0;
    tmpParams.paramMainMatrix([5,6])=0;
    tmpParams.paramMainMatrix([4])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .25 .25],'NewWindow',false,'DisplayName', 'transient');
    pbaspect([1 2 1])
    xlim([0 13]);
    set(h,'Xtick',0:2:13)
    hold on
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,6])=0;
    tmpParams.paramMainMatrix([5])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .5 .5],'NewWindow',false,'DisplayName', 'sustained');
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,5])=0;
    tmpParams.paramMainMatrix([6])=1;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .75 .75],'NewWindow',false,'DisplayName', 'persistent');
    legend('show')
    hold off
    % Plot the scaled final components
    h=subplot(2,2,4);
    tmpParams=params0;
    tmpParams.paramMainMatrix([5,6])=0;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .25 .25],'NewWindow',false,'DisplayName', 'transient');
    pbaspect([1 2 1])
    xlim([0 13]);
    set(h,'Xtick',0:2:13)
    hold on
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,6])=0;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .5 .5],'NewWindow',false,'DisplayName', 'sustained');
    tmpParams=params0;
    tmpParams.paramMainMatrix([4,5])=0;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[1 .75 .75],'NewWindow',false,'DisplayName', 'persistent');
    tmpParams=params0;
    modelResponseStruct=temporalFit.computeResponse(tmpParams,stimulusStruct,[]);
    temporalFit.plot(modelResponseStruct,'Color',[.75 .75 .75],'NewWindow',false,'DisplayName', 'fullm model');
    legend('show')
    hold off
end

%% Set returned validationData structure
if (nargout > 0)
    validationData.params1 = paramsFit;
    validationData.modelResponseStruct = modelResponseStruct;
end

end
