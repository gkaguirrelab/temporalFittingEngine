function packetsConc = concatenatePackets(obj,thePackets)

% function packetsConc = concatenatePackets(packets)
%
% takes in a list of packets, and concatenates their stimulus and response
% timebases, as well as their stimulus and response values. returns a
% struct.

p = inputParser;
p.addRequired('thePackets',@iscell);
p.parse(thePackets);

% initialize struct to return
packetsConc = struct;

% vector for concatenating time series'
responseValues = [];
% vector for concatenating stimulus matrices
stimulusValues = [];
% response and stimulus time bases
responseTB = [];
stimulusTB = [];

% determine the number of events in the packet with the most events
numEventsStore = [];
% determine the total length of the run
lengthRun = 0;
for i = 1:length(thePackets)
   numEventsStore(i) = size(thePackets{i}.stimulus.values,1);
   lengthRun = lengthRun + size(thePackets{i}.stimulus.values,2);
end
numEvents = max(numEventsStore);
%%

% mark the index of the stimulus values
packetPositionMarker = 1;
% mark the end of the current stimulus and response timebases
stimulusTBmarker = 0;
responseTBmarker = 0;

for i = 1:length(thePackets)
    % concatenate response values
    responseValues = [responseValues thePackets{i}.response.values];
    % get the current packet's stimulus values. start by creating a matrix
    % of 0's that has the number of stimuli in the packet as the first
    % dimension, and the total length of the concatenated stimulus as the
    % second
    stimulusValuesForOnePacket = zeros([size(thePackets{i}.stimulus.values,1) lengthRun]);
    % stick the current packet's stimulus values in this matrix of zeroes
    stimulusValuesForOnePacket(1:size(thePackets{i}.stimulus.values,1), ...
                               packetPositionMarker:packetPositionMarker+ ...
                               size(thePackets{i}.stimulus.values,2)-1) ...
                               = thePackets{i}.stimulus.values;
    % update the index marker of the packet stimulus values
    packetPositionMarker = packetPositionMarker+size(thePackets{i}.stimulus.values,2);
    % concatenate stimulus values
    stimulusValues = [stimulusValues; stimulusValuesForOnePacket];
    
    % response timebase for one packet: time-shift by timebase marker
    responseTBforOnePacket = thePackets{i}.response.timebase+responseTBmarker;
    % concatenate
    responseTB = [responseTB responseTBforOnePacket];
    
    % do the same for stimulus
    stimulusTBforOnePacket = thePackets{i}.stimulus.timebase+stimulusTBmarker;
    stimulusTB = [stimulusTB stimulusTBforOnePacket];
    
    % update timebase markers
    stimulusTBmarker = stimulusTB(length(stimulusTB));
    responseTBmarker = responseTB(length(responseTB));
    
end

stimulusTB = linspace(min(stimulusTB),max(stimulusTB),length(stimulusTB));
responseTB = linspace(min(responseTB),max(responseTB),length(responseTB));

% assign to struct
packetsConc.stimulus.timebase = stimulusTB;
packetsConc.stimulus.values = stimulusValues;
packetsConc.response.timebase = responseTB;
packetsConc.response.values = responseValues;

end