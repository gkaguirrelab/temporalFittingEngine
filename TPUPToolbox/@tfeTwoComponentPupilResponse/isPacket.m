function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.
%
% 9/14/16   ms      Wrote it.

if isfield(thePacket, 'stimulus') && ...
        isfield(thePacket, 'response') && ...
        isfield(thePacket, 'kernel') && ...
        isfield(thePacket, 'metaData') && ...
        length(thePacket) == 1
    packetValidity = true;
else
    packetValidity = false;
end