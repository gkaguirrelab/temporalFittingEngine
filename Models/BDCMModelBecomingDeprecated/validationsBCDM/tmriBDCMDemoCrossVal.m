%% LOAD PACKETS

load('/Users/benchin/Desktop/MELA_data/packets_asb1_041416.mat');
%% TRY CROSS-VALIDATION FOR JUST THE L-M
indicesOfInterest = [2 5 8 11];

thePackets = {};

for i = 1:length(indicesOfInterest)
   aPacket = packets{indicesOfInterest(i)};
   aPacket.response.values = aPacket.response.values{indicesOfInterest(i)};
   thePackets{length(thePackets)+1} = aPacket; 
end

clear packets;

%% Construct the model object
tmri = tmriBDCM;

%% DOWNSAMPLE STIMULI TO THE DESIRED RESOLUTION

% common factor to downsample by
timeFactor = 200;

startTime = thePackets{1}.stimulus.timebase(1);
endTime = thePackets{1}.stimulus.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newStimulusTimebase = linspace(startTime,endTime,nTimeSamples+1);

%% DOWNSAMPLE HRF TO THE DESIRED RESOLUTION

startTime = thePackets{1}.HRF.timebase(1);
endTime = thePackets{1}.HRF.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newHRFTimebase = linspace(startTime,endTime,nTimeSamples+1);

%%
% May want a utility function to convert lists of packets into lists of
% their substructures, and back.
for ii = 1:length(thePackets)
    theStimulusList{ii} = thePackets{ii}.stimulus;
    theHRFList{ii} = thePackets{ii}.HRF;
end

resampledStimulusList = tmri.resamplePacketStruct(theStimulusList,newStimulusTimebase);
resampledHRFList = tmri.resamplePacketStruct(theHRFList,newHRFTimebase);

for ii = 1:length(thePackets)
    thePackets{ii}.stimulus = resampledStimulusList{ii};
    theHRFList{ii} = resampledHRFList{ii};
    thePackets{ii}.HRF = theHRFList{ii};
end

%% CROSS VALIDATION, FITTING ALL BUT PACKET #3

[fVal responsePredicted] = tmri.crossValBDCM(thePackets,3);