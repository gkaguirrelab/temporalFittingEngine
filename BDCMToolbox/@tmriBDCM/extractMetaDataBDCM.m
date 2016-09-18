function metaData = extractMetaDataBDCM(obj,packet)

% function metaData = extractMetaDataBDCM(obj,packet)
%
% for BDCM study. Takes in packet, and pulls out frequency indices

p = inputParser;
p.addRequired('packet',@isstruct);
p.parse(packet);

metaData = struct;

metaData.theFrequencyIndices = packet.stimulus.metaData.params.theFrequencyIndices;

end