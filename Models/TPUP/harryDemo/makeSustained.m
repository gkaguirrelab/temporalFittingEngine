function makeSustained(exponentialTauVec)

gammaTauVec = 200;

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

%% now get the forward model going
obj = tfeTPUP('verbosity','full');

gammaIRF.values = stimulusStruct.timebase .* 0;
gammaIRF.timebase = stimulusStruct.timebase;
exponentialIRF.values=stimulusStruct.timebase .* 0;
exponentialIRF.timebase=stimulusStruct.timebase;

ii = 1;
stimulus.values = stimulusStruct.values(ii,:);
stimulus.timebase = stimulusStruct.timebase;

gammaIRF.values = stimulus.timebase .* exp(-stimulus.timebase./gammaTauVec(ii));
gammaIRF=normalizeKernelArea(gammaIRF);


exponentialIRF.values=exp(-1/exponentialTauVec(ii)*stimulus.timebase);
exponentialIRF=normalizeKernelArea(exponentialIRF);

sustainedComponent = obj.applyKernel(stimulus,gammaIRF);
sustainedComponent=normalizeKernelArea(sustainedComponent);


convolvedSustainedComponent = obj.applyKernel(obj.applyKernel(stimulus,exponentialIRF),gammaIRF);
convolvedSustainedComponent=normalizeKernelArea(convolvedSustainedComponent);

plot(sustainedComponent.timebase, convolvedSustainedComponent.values)
