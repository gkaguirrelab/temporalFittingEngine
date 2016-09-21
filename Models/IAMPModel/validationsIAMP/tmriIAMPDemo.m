% tmriInstanceAmplitudeDemo
%
% Demonstrate the operation of a very simple temporal fitting model
%
% 09/13/16  gka  Wrote it.

%% Clear and close
clear; close all;

%% Construct the model object
temporalFit = tmriIAMP;

%% Load in a packet to fit
%
% The packets variable in the .mat file is a cell array of packets, each a
% structure.
theExampleData = load('/Users/aguirre/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/HCLV_photo_7T/packets_asb1_041416.mat');
nPacketsRead = length(theExampleData.packets);

% For the demo, just take the first packet and simplify the variable name
nPackets = 1;
thePackets = cell(nPackets,1);
for ii = 1:length(thePackets)
    thePackets{ii} = theExampleData.packets{ii};
end
clear theExampleData

%Fix example packets to look like what we have decided real packets will
%be like.
for ii = 1:length(thePackets)
    thePackets{ii}.stimulus.timebase = thePackets{ii}.stimulus.timebase(1,:);
    thePackets{ii}.response.timebase = thePackets{ii}.response.timebase(1,:);
    thePackets{ii}.HRF.timebase = thePackets{ii}.HRF.timebase(1,:);
    
%     for ct = 1:length(thePackets{ii}.metaData.fileName)
% %        thePackets{ii}.stimulus.metaData(ct).fileName = thePackets{ii}.metaData.fileName{ct};
%         thePackets{ii}.stimulus.metaData(ct).frequency = 1; %BDCMfilenameReader(thePackets{ii}.metaData.fileName{ct});
%         thePackets{ii}.stimulus.metaData(ct).colordir = 'Foo'; % Figure me out
%     end
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
resampledStimulusList = temporalFit.resamplePacketStruct(theStimulusList,newStimulusTimebase);
resampledHRFList = temporalFit.resamplePacketStruct(theHRFList,newHRFTimebase);
for ii = 1:nPackets
    thePackets{ii}.stimulus = resampledStimulusList{ii};
    theHRFList{ii} = resampledHRFList{ii};
    thePackets{ii}.HRF = theHRFList{ii};
end

%% How many stimulus instances were there
defaultParamsInfo.nEvents = size(thePackets{1}.stimulus.values,1);

%% Obtain default parameters
%
params0 = temporalFit.defaultParams('DefaultParamsInfo',defaultParamsInfo);
fprintf('Default model parameters:\n');
temporalFit.print(params0);

%% Test paramsToVec and vecToParams
params1 = params0;
params1.paramsMainMatrix = rand(size(params1.paramMainMatrix));
x1 = temporalFit.paramsToVec(params1);
params2 = temporalFit.vecToParams(x1);
if (any(params1.paramMainMatrix ~= params2.paramMainMatrix))
    error('vecToParams and paramsToVec do not invert');
end

%% Set a parameter locking matrix
% The empty matrix is set here, so no locking
paramLockMatrix = [];

%% Conduct the fit
[paramsFit,fVal,fitResponse] = temporalFit.fitResponse(thePackets{1},'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
fprintf('Model parameter from fits:\n');
oneResponse = fitResponse{1};

