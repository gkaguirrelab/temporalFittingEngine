function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

% Basic structure of a test. Replace the test here with something
% appropriate for this model.
if ~numel((thePacket.stimulus.values)==3)
    warning('The pRF model requires a movie stimulus (3 dimensions)')
    packetValidity = false;
end