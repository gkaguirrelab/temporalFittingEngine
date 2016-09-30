function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.

% the packet is innocent until proven guilty
packetValidity = true;

% let the user know what we are about to do
switch (obj.verbosity)
    case 'high'
        fprintf('Checking the passed packet...');
end

% Check for the presence of primary fields
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

% Check for the presence of secondary fields
% Both the stimulus and response fields must have a timebase and values.
% metaData is optional at the secondary level. Definition of the kernel is
% optional.
if isfield(thePacket.stimulus, 'values') && ...
        isfield(thePacket.stimulus, 'timebase') && ...
        isfield(thePacket.response, 'values') && ...
        isfield(thePacket.response, 'timebase')
    packetValidity = true;
else
    warning('Either timebase or values are missing from stimulus or response')
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
if ~(size(thePacket.stimulus.timebase,1)==1)
    warning('The field response.timebase is not a single row vector')
    packetValidity = false;
end

% Exit at this stage if the packet is bad, as we can't check the
% response field lengths
if ~packetValidity
    return
end

% Check that values and timebase are the same length for the response
if ~(length(thePacket.response.values)==length(thePacket.response.timebase))
    warning('response.timebase is not equal in length to response.values')
    packetValidity = false;
end

% Check that the second dimension of values and timebase are the same
%  length for the stimulus
if ~(size(thePacket.stimulus.values,2)==length(thePacket.stimulus.timebase))
    warning('response.timebase is not equal in length to response.values')
    packetValidity = false;
end

% If the kernel values and timebase fields are defined and not empty,
%  make sure that they are of the proper dimensions and of the same length
if isfield(thePacket.kernel, 'values') && ...
        isfield(thePacket.kernel, 'timebase')
    if ~isempty(thePacket.kernel.values) &&...
            ~isempty(thePacket.kernel.timebase)

        % the kernel fields are defined and not empty. Check them.
        if ~(size(thePacket.kernel.values,1)==1)
            warning('The field kernel.value is defined but not a single row vector')
            packetValidity = false;
        end
        if ~(size(thePacket.kernel.timebase,1)==1)
            warning('The field kernel.timebase is defined but not a single row vector')
            packetValidity = false;
        end
        if ~(size(thePacket.kernel.values)==length(thePacket.kernel.timebase))
            warning('kernel.timebase and kernel.values are defined but not equal in length')
            packetValidity = false;
        end
        
    end % if the kernel fields are not empty
end % if the kernel fields are defined

% Check that the stimulus fields are not empty
if isempty(thePacket.stimulus.values)
    warning('The field stimulus.values is empty')
    packetValidity = false;
end
if isempty(thePacket.stimulus.timebase)
    warning('The field stimulus.timebase is empty')
    packetValidity = false;
end

% Check that the response fields are not empty
if isempty(thePacket.response.values)
    warning('The field response.values is empty')
    packetValidity = false;
end
if isempty(thePacket.response.timebase)
    warning('The field response.timebase is empty')
    packetValidity = false;
end
