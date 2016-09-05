% tmriTFBlockDesignColorModelDemo
%
% Demonstrate function for the quadratic color model.
%
% 6/26/16  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Construct the model object
tmri = tmriTFBlockDesignColorModel;

%% Load in some packets to play with
%
% Eventually, we will update this to use a fancier packet delivery
% method, but for now we just want to prepare for that joyous day.
%
% The packets variable in the .mat file is a cell array of packets, each a
% structure.
theExampleData = load(fullfile('BDCMTestData','packets_asb1_041416.mat'));
nPacketsRead = length(theExampleData.packets);
nPackets = 1;
thePackets = cell(nPackets,1);
for ii = 1:length(thePackets)
    thePackets{ii} = theExampleData.packets{ii};
end
clear theExampleData

% Fix example packets to look like what we have decided real packets will
% be like.
for ii = 1:length(thePackets)
    thePackets{ii}.stimulus.timebase = thePackets{ii}.stimulus.timebase(1,:);
    thePackets{ii}.response.timebase = thePackets{ii}.response.timebase(1,:);
    thePackets{ii}.HRF.timebase = thePackets{ii}.HRF.timebase(1,:);
    
    for ct = 1:length(thePackets{ii}.metaData.fileName)
        thePackets{ii}.stimulus.metaData(ct).fileName = thePackets{ii}.metaData.fileName{ct};
        thePackets{ii}.stimulus.metaData(ct).frequency = BDCMfilenameReader(thePackets{ii}.metaData.fileName{ct});
        thePackets{ii}.stimulus.metaData(ct).colordir = 'Foo'; % Figure me out
    end
end

%% Downsample timebase to something more manageable
%
% Time is in milliseconds, so with timeFactor = 100 we go to 
% roughly 100 millisecond time steps.
%
% Here we are assuming that all packet stimuli live on the same
% timebase as each other, and same for all packet HRFs.  Probably should
% check or generalize.
timeFactor = 100;
startTime = thePackets{1}.stimulus.timebase(1);
endTime = thePackets{1}.stimulus.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newStimulusTimebase = linspace(startTime,endTime,nTimeSamples+1);
startTime = thePackets{1}.HRF.timebase(1);
endTime = thePackets{1}.HRF.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newHRFTimebase = linspace(startTime,endTime,nTimeSamples+1);

% May want a utility function to convert lists of packets into lists of
% their substructures, and back.
for ii = 1:nPackets
    theStimulusList{ii} = thePackets{ii}.stimulus;
    theHRFList{ii} = thePackets{ii}.HRF;
end
resampledStimulusList = tmri.resamplePacketStruct(theStimulusList,newStimulusTimebase);
resampledHRFList = tmri.resamplePacketStruct(theHRFList,newHRFTimebase);
for ii = 1:nPackets
    thePackets{ii}.stimulus = resampledStimulusList{ii};
    theHRFList{ii} = resampledHRFList{ii};
    thePackets{ii}.HRF = theHRFList{ii};
end

%% How many individual stimuli were there
%
% This needs a little thought and perhaps should become a makeDesign sort
% of function
defaultParamsInfo.nEvents = size(thePackets{1}.stimulus.values,1);

%% Specify the stimulus and response. 
%
% The information is currently spread across two test structures, so we
% load them both an re-arrange to our liking here.
%
% Be sure to define defaultParamsInfo.nEvents here, as we need it to create default
% parameters for the model.
% whichRunToTest = 1;
% stimulus.stimMatrix = squeeze(theTestStimuli.stimulusStruct.stimMatrix(whichRunToTest,:,:));
% stimulus.stimValues = squeeze(theTestStimuli.stimulusStruct.stimValuesForRunStore(whichRunToTest,:));
% stimulus.startTimesSorted_A = theTestStimuli.stimulusStruct.startTimesSorted_A;
% stimulus.startTimesSorted_B = theTestStimuli.stimulusStruct.startTimesSorted_B;
% stimulus.stimValuesSorted_A = theTestStimuli.stimulusStruct.stimValuesSorted_A;
% stimulus.stimValuesSorted_B = theTestStimuli.stimulusStruct.stimValuesSorted_B;
% stimulus.uniqueTemporalFreq = theTestStimuli.stimulusStruct.uniqueTemporalFreq;
% stimulus.AOrB = theTestData.exampleResponses.runOrder(whichRunToTest);
% stimulus.modulationDir = theTestData.exampleResponses.modulationDir(whichRunToTest);
% boldResponse = squeeze(theTestData.exampleResponses.cleanedData(whichRunToTest,:));
% clearvars('theTestStimuli','theTestData');

%% Set parameters
%
% Six parameters define a quadratic form in three dimensions, but
% we normalize the first to 1 so we only need five numbers here.
params0 = tmri.defaultParams('DefaultParamsInfo',defaultParamsInfo);
fprintf('Default model parameters:\n');
tmri.print(params0);

%% Test paramsToVec and vecToParams
params1 = params0;
params1.paramsMainMatrix = rand(size(params1.paramMainMatrix));
x1 = tmri.paramsToVec(params1);
params2 = tmri.vecToParams(x1);
if (any(params1.paramMainMatrix ~= params2.paramMainMatrix))
    error('vecToParams and paramsToVec do not invert');
end

%% Get parameter locking matrix
% paramLockMatrix = tmri.lockMatrix(params0,thePacket.stimulus);
paramLockMatrix = [];

%% Plot the BOLD response
%
% Need to modify this to handle lists of packets.  Should we pass the lists
% to plot, or should we loop around the plot function?
% tmri.plot(thePacket.response.timebase,thePacket.response.values);

%% Test the fitter
% [paramsFit,fitResponse] = tmri.fitResponse(timebase,stimulus,boldResponse, ...
%                           'HRF',theHRF,'DefaultParamsInfo',defaultParamsInfo);
[paramsFit,fVal,fitResponse] = tmri.fitResponse(thePackets,'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
fprintf('Model parameter from fits:\n');
oneResponse = fitResponse{1};
tmri.plot(thePackets{1}.response.timebase,thePackets{1}.response.values); hold on;
tmri.plot(thePackets{1}.response.timebase,oneResponse{1},'Color',[0 1 0],'NewWindow',false);
%tmri.print(paramsFit);
%tmri.plot(thePacket.stimulus.timebase(1,:),fitResponse{1},'Color',[0 1 0],'NewWindow',false);
% [~,meanParamValues] = tmri.plotParams(paramsFit,stimulus);

%% Test that we can obtain a neural response
% fprintf('Simulated model parameters:\n');
% tmri.print(params1);
% responseToFit = tmri.computeResponse(params1,timebase,stimulus,'AddNoise',true);
