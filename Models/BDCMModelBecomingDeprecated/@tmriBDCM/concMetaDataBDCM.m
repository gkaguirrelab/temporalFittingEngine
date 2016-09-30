function metaData = concMetaDataBDCM(obj,thePackets)

% function metaData = concMetaDataBDCM(obj,thePackets)
% 
% loops through a bunch of packets and extracts metadata for BDCM, then
% concatenates it. Returns a struct.

p = inputParser;
p.addRequired('thePackets',@iscell);
p.parse(thePackets);

theFrequencyIndicesStore = [];

tmri = tmriBDCM;

for i = 1:length(thePackets)
   curMetaData = tmri.extractMetaDataBDCM(thePackets{i});
   theFrequencyIndicesStore = [theFrequencyIndicesStore curMetaData.theFrequencyIndices];
end

metaData.theFrequencyIndices = theFrequencyIndicesStore;

end