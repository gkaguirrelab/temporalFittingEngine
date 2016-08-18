% tmriTFBlockDesignColorModelDemo
%
% Demonstrate function for the quadratic color model.
%
% 6/26/16  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'..'));
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'toolbox'));

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
nPackets = 4;
thePackets = theExampleData.packets(1:nPackets);
clear theExampleData

%% Start by analyzing just one packet
%
% Adjust the example packet format to match the format that will eventually
% be supplied
thePacket = thePackets{1};
for ct = 1:length(thePacket.metaData.fileName)
    thePacket.stimulus.metaData(ct).fileName = thePacket.metaData.fileName{ct};
    thePacket.stimulus.metaData(ct).frequency = 0; % Figure me out
    thePacket.stimulus.metaData(ct).colordir = 'Foo'; % Figure me out
end

%% How many individual stimuli were there
defaultParamsInfo.nEvents = size(thePacket.stimulus.values,1);

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

%% Set the timebase we want to compute on
%
% For now, we know that the TRs are in one second units, and we have
% decided to work in one second units.
% deltaT = 1;
% totalTime = size(stimulus.stimMatrix,2)-1;
% timebase = 0:deltaT:totalTime;

%% Get the HRF
% theTestHRF = load(fullfile('BDCMTestData','asb1HRF.mat'));
% theHRF.HRF = theTestHRF.BOLDHRF;
% theHRF.timebase = 0:deltaT:(size(theHRF.HRF,2)-1);
% clearvars('theTestHRF');

%% Get parameter locking matrix
% paramLockMatrix = tmri.lockMatrix(params0,stimulus);
paramLockMatrix = [];

%% Plot the BOLD response
tmri.plot(thePacket.response.timebase,thePacket.response.values);

%% Test the fitter
% [paramsFit,fitResponse] = tmri.fitResponse(timebase,stimulus,boldResponse, ...
%                           'HRF',theHRF,'DefaultParamsInfo',defaultParamsInfo);
[paramsFit,fVal,allFVals,fitResponse] = tmri.fitResponse({thePacket.stimulus},{thePacket.response}, ...
                          'HRF',{thePacket.HRF},'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
fprintf('Model parameter from fits:\n');
tmri.print(paramsFit);
tmri.plot(thePacket.stimulus.timebase(1,:),fitResponse{1},'Color',[0 1 0],'NewWindow',false);
% [~,meanParamValues] = tmri.plotParams(paramsFit,stimulus);

%% Test that we can obtain a neural response
% fprintf('Simulated model parameters:\n');
% tmri.print(params1);
% responseToFit = tmri.computeResponse(params1,timebase,stimulus,'AddNoise',true);
