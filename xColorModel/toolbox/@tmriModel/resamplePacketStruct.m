function resampledStructList = resamplePacketStruct(obj,packetStructList,newTimebase,varargin)
% function resampledStruct = resamplePacketStruct(packetStructList,newTimebase,varargin
%
% Resamples timebase and values within a list of structs to the passed
% new timebase.  Preserves the other fields of each struct in the passed
% list.
%
% Required Inputs:
%   packetStructList: Cell array of structs with timebase and values fields.
% 	  Typically these will be elemements of a packet.
%   newTimebase: New timebase
%
% Optional key/value pairs:
%   'Method'
%     'interp1_linear' (default) - Use Matlab's interp1 linear method.

% Parse input
p = inputParser;
p.addRequired('packetStructList',@iscell);
p.addRequired('newTimebase',@isnumeric);
p.addParameter('Method','interp1_linear',@ischar);
p.parse(packetStructList,newTimebase,varargin{:});

% Loop over the structs in the cell array.
%
% For each, duplicate the input structure, but empty the timebase and values fields
nStructs = length(packetStructList);
resampledStructList = cell(size(packetStructList));
for ll = 1:nStructs
    resampledStruct = packetStructList{ll};
    resampledStruct.timebase = newTimebase;
    resampledStruct.values = [];
    
    % Loop over each row of values, and downsample
    for ii = 1:size(packetStructList{ll}.values,1)
        switch (p.Results.Method)
            case 'interp1_linear'
            	resampledStruct.values(ii,:) = interp1(packetStructList{ll}.timebase,packetStructList{ll}.values(ii,:),newTimebase,'linear');
            otherwise
                error('Unknown method specified');
        end
    end
    
    % Tuck this into the retun list
    resampledStructList{ll} = resampledStruct;
end
end