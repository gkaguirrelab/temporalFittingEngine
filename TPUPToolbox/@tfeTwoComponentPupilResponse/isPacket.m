function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.
%
% 9/14/16   ms      Wrote it.
% 9/21/16   gka     Expanded functionality and reporting

% the packet is innocent until proven guilty
packetValidity = true;

% let the user know what we are about to do
fprintf('Checking the passed packet...\n');

% Check for the presence of elementary fields
if isfield(thePacket, 'stimulus') && ...
        isfield(thePacket, 'response') && ...
        isfield(thePacket, 'kernel') && ...
        isfield(thePacket, 'metaData') && ...
        length(thePacket) == 1
    packetValidity = true;
else
    warning('An elementary packet field is missing')
    packetValidity = false;
end

% exit at this stage if the packet is bad, as we can't check the fields
if ~packetValidity
    return
end

% Check for proper dimensionality of the response fields
if ~(size(thePacket.response.values,1)==1)
    warning('The field response.value is not a single row vector')
    packetValidity = false;
end
if ~(size(thePacket.response.timebase,1)==1)
    warning('The field response.timebase is not a single row vector')
    packetValidity = false;
end

% exit at this stage if the packet is bad, as we can't check the the
% response field lengths
if ~packetValidity
    return
end

% Check that the response values and timebase are the same length
if ~(length(thePacket.response.values)==length(thePacket.response.timebase))
    warning('response.timebase is not equal in length to response.values')
    packetValidity = false;
end

% Check that the stimulus fields are not empty
if isempty(thePacket.stimulus.values)
    warning('The field stimulus.values is empty')
    packetValidity = false;
end
if isempty(thePacket.stimulus.timebase)
    warning('The field stimulus.timebase is empty')
    packetValidity = false;
end

% Check to see if the initial value of response.values is a small number
if thePacket.response.values(1) > 0.001
    warning('The first value of response.values is not approximately zero')
    packetValidity = false;
end
