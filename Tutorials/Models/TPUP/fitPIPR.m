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

%% Grab net PIPR result
[ BlueAverage ] = averageResultAcrossSubjects(averageBlueCombined{1});
[ RedAverage ] = averageResultAcrossSubjects(averageRedCombined{1});

thePacket.response.values = BlueAverage - RedAverage;
thePacket.response.timebase = timebase;

%% now actually do the fit
initialValues=[200, 200, 12000, -10, -25, -25];
[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(thePacket, 'defaultParamsInfo', defaultParamsInfo, 'initialValues', initialValues);

figure;
hold on
plot(timebase, BlueAverage, 'Color', 'b')
plot(timebase, RedAverage, 'Color', 'r')
plot(timebase, thePacket.response.values, 'Color', 'k')
plot(timebase, modelResponseStruct.values, 'Color', 'c')
xlabel('Time (ms)')
ylabel('Pupil Diameter (% Change)')
legend('Blue Response Average', 'Red Response Average', 'PIPR Average', 'Model Fit', 'Location', 'SouthEast')


%plot fit with each component
% 1st component
paramsFitTemp = paramsFit;
paramsFitTemp.paramMainMatrix(5) = 0;
paramsFitTemp.paramMainMatrix(6) = 0;
tmp1 = temporalFit.computeResponse(paramsFitTemp,stimulusStruct,thePacket.kernel,'AddNoise',false);

% 2nd component
paramsFitTemp = paramsFit;
paramsFitTemp.paramMainMatrix(4) = 0;
paramsFitTemp.paramMainMatrix(6) = 0;
tmp2 = temporalFit.computeResponse(paramsFitTemp,stimulusStruct,thePacket.kernel,'AddNoise',false);

% 3rd component
paramsFitTemp = paramsFit;
paramsFitTemp.paramMainMatrix(4) = 0;
paramsFitTemp.paramMainMatrix(5) = 0;
tmp3 = temporalFit.computeResponse(paramsFitTemp,stimulusStruct,thePacket.kernel,'AddNoise',false);

figure;
hold on
plot(timebase, tmp1.values)
plot(timebase, tmp2.values)
plot(timebase, tmp3.values)
xlabel('Time (ms)')
ylabel('Component Amplitude')
legend('Transient', 'Sustained', 'Persistent', 'Location', 'SouthEast')
