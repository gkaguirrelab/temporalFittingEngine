function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.
%
% 9/14/16   ms      Wrote it.
% 9/21/16   gka     Expanded functionality and reporting

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

% Check to see if the initial value of response.values is a small number
if thePacket.response.values(1) > 0.001
    warning('The first value of response.values is not approximately zero')
    packetValidity = false;
end
