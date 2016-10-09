function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

% Basic structure of a test. Replace the test here with something
% appropriate for this model.
if ~packetValidity
    warning('Yup. Still bad.')
    packetValidity = false;
end
