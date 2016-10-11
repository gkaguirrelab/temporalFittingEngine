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
        isfield(thePacket.stimulus.fcon, 'observedParamMatrix') && ...
        isfield(thePacket.stimulus.fcon, 'modelObjHandle') && ...
        isfield(thePacket.stimulus.fcon, 'logContrastFlag')
    packetValidity = true;
else
    warning('There are fields missing from packet.stimulus.fcon')
    packetValidity = false;
    return
end

% Test if the dimension of contrastbase matches that of observedParamMatrix
if ~(length(thePacket.stimulus.fcon.contrastbase)==size(thePacket.stimulus.fcon.observedParamMatrix,2))
    warning('the contrastbase must be the same length as the second dimension of the observedParamMatrix')
    packetValidity = false;
end

% Test if the logContrastFlag is boolean
if ~isboolean(thePacket.stimulus.fcon.logContrastFlag)
    warning('the logContrastFlag field does not contain a boolean value')
    packetValidity = false;
end

% Test if the modelObjHandle is an object handle (Not sure how to implement
% yet)
