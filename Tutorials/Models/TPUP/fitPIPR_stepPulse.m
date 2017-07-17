%% fitPIPR

%% basics
defaultParamsInfo.nInstances = 1;

% Construct the model object
temporalFit = tfeTPUP('verbosity','full');

%% make stimulus structure
% assemble the packet
% first create the stimulus structure

% create the timebase: events are 14 s long, and we're sampling every 20
% ms
timebase = (0:20:13998);




% Temporal domain of the stimulus
deltaT = 20; % in msecs
totalTime = 14000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

% Specify the stimulus struct.
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
thePacket.stimulus.values = stimulusStruct.values;
thePacket.stimulus.timebase = timebase;

% now kernel needed for tpup
thePacket.kernel = [];
thePacket.metaData = [];

%% first play around with making a dummy fit
thePacket.response.values = zeros(1,length(timebase));
thePacket.response.values(50:600) = -10;
thePacket.response.timebase = timebase;

plotFig = figure;
hold on
plot(thePacket.response.timebase, thePacket.response.values)

[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(thePacket, 'defaultParamsInfo', defaultParamsInfo);

plot(modelResponseStruct.timebase, modelResponseStruct.values)
paramsFit
