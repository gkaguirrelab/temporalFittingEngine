function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

% Check to see if the initial value of response.values is a small number.
% This is instantiated by checking if the first value is less than |1%| of
% the entire range of the values.
packetValuesRange=max(thePacket.response.values)-min(thePacket.response.values);
if abs(thePacket.response.values(1) / packetValuesRange) > 0.01
    warning('The first value of response.values is not approximately zero')
    packetValidity = false;
end
