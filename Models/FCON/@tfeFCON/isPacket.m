function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

% Test for the presence of the fcon fields within the stimulusStruct portion
% of the packet
if isfield(thePacket.stimulus, 'fcon')
    packetValidity = true;
else
    warning('The FCON model requires a the field packet.stimulus.fcon')
    packetValidity = false;
    return
end

if isfield(thePacket.stimulus.fcon, 'contrastbase') && ...
        isfield(thePacket.stimulus.fcon, 'paramLookUpMatrix') && ...
        isfield(thePacket.stimulus.fcon, 'modelObjHandle') && ...
        isfield(thePacket.stimulus.fcon, 'defaultParams')
    packetValidity = true;
else
    warning('There are fields missing from packet.stimulus.fcon')
    packetValidity = false;
    return
end

% Test if the dimension of contrastbase matches that of paramLookUpMatrix
if ~(length(thePacket.stimulus.fcon.contrastbase)==size(thePacket.stimulus.fcon.paramLookUpMatrix,2))
    warning('The contrastbase must be the same length as the second dimension of the paramLookUpMatrix')
    packetValidity = false;
end

% Test that the defaultParams are a structure
if ~(isstruct(thePacket.stimulus.fcon.defaultParams))
    warning('The defaultParams field does not contain a struct')
    packetValidity = false;
end

% Test if the modelObjHandle is an object handle (Not sure how to implement
% yet)
