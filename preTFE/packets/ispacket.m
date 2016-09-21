function out = ispacket(packet)
% out = ispacket(packet)
%
% Function to test if the passed `packet` is well-formed.
%
% 9/14/16   ms      Wrote it.

if isfield(packet, 'stimulus') && ...
        isfield(packet, 'response') && ...
        isfield(packet, 'kernel') && ...
        isfield(packet, 'metaData') && ...
        length(packet) == 1
    out = true;
else
    out = false;
end