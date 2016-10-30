function [ theConcatPacket ] = concatenatePackets(obj,packetCellArray, varargin)
% function theConcatPacket = concatenatePackets(obj,packetCellArray)
%
% Concatenates a cell array of packets.
%
% If stimLabels are defined in the first packet, then all packets are
% required to have stimLabels and stimTypes. It is further then required
% that all packets have the same set of stimLabels, in the same order. The
% different packets may have stimTypes that do not include instances of all
% of the available stimLabels.
%
% If defined, all the packets must have identical kernel fields.
%
%

%% Parse vargin for options passed here
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('packetCellArray',@iscell);
p.addParameter('stimValueExtender',0,@isnumeric);
p.parse(packetCellArray, varargin{:});

%% Check the input
% Determine the number of packets
nPackets=length(packetCellArray);

% If there is only one packet, set theConcatPacket to this and return
if nPackets == 1
    theConcatPacket=packetCellArray;
    return
end

% If the kernel field exists, make sure it is identical across the packets
if isfield(packetCellArray{1},'kernel')
    holdKernelStruct=packetCellArray{1}.kernel;
    for pp=2:nPackets
        if ~isequal(holdKernelStruct,packetCellArray{pp}.kernel)
            error('All packets must contain the identical kernel struct');
        end % check for identical kernel Structs
    end % loop over packets
end

% if the first packet has the field stimLabels, then check that every
% packet has the same stimulus labels, and that all packets have stimTypes
if isfield(packetCellArray{1}.stimulus.metaData,'stimLabels')
    if ~isfield(packetCellArray{1}.stimulus.metaData,'stimTypes')
        error('If stimLabels are defined, all packets must have stimTypes');
    end % check the first packet for stimTypes
    uniqueStimLabels=unique(packetCellArray{1}.stimulus.metaData.stimLabels);
    for pp=2:nPackets
        if ~isfield(packetCellArray{pp}.stimulus.metaData,'stimTypes')
            error('If stimLabels are defined, all packets must have stimTypes');
        end % check the first packet for stimTypes
        if  ~isfield(packetCellArray{pp}.stimulus.metaData,'stimLabels')
            error('Not all of the packets have a stimLabels field');
        end % test for the existence of stimLabels field
        if  ~isequal(uniqueStimLabels, unique(packetCellArray{pp}.stimulus.metaData.stimLabels))
            error('The packets do not have identical stimLabels');
        end % test for the same stimLabels
    end
end % the first packet has a stimLabel


%% Build the concatenated packet

% Initialize fields for theConcatPacket
theConcatPacket.stimulus.timebase=[];
theConcatPacket.stimulus.metaData.stimTypes=[];
theConcatPacket.response.timebase=[];
theConcatPacket.response.values=[];

% If defined, copy over the stimLabels
if isfield(packetCellArray{1}.stimulus.metaData,'stimLabels')
    theConcatPacket.stimulus.metaData.stimLabels= ...
        packetCellArray{1}.stimulus.metaData.stimLabels;
end

% Loop throught the packets and assemble the concatenated timebase
stimTimeMarker=0;
respTimeMarker=0;
totalStimValuesRows=0;
for pp=1:nPackets
    
    % Add up the total number of stim value rows. We'll need this later
    totalStimValuesRows=totalStimValuesRows+size(packetCellArray{pp}.stimulus.values,1);
    
    % Concatenate the response vector
    theConcatPacket.response.values = [theConcatPacket.response.values ...
        packetCellArray{pp}.response.values];
    
    % If defined, concatenate the stimTypes
    if isfield(packetCellArray{1}.stimulus.metaData,'stimTypes')
        theConcatPacket.stimulus.metaData.stimTypes = ...
            [theConcatPacket.stimulus.metaData.stimTypes; ...
            packetCellArray{pp}.stimulus.metaData.stimTypes];
    end
    
    % Concatenate the timebase, advancing time as we go. Assume that the
    % time elapsed between the end of one packet and the start of the next
    % is the same as the time between the last two timebase values
    theConcatPacket.stimulus.timebase = ...
        [theConcatPacket.stimulus.timebase ...
        packetCellArray{pp}.stimulus.timebase+stimTimeMarker];
    stimTimeMarker=theConcatPacket.stimulus.timebase(end) + ...
        (theConcatPacket.stimulus.timebase(end) - theConcatPacket.stimulus.timebase(end-1));
    theConcatPacket.response.timebase = ...
        [theConcatPacket.response.timebase ...
        packetCellArray{pp}.response.timebase+respTimeMarker];
    respTimeMarker=theConcatPacket.response.timebase(end) + ...
        (theConcatPacket.response.timebase(end) - theConcatPacket.response.timebase(end-1));
    
    % place all other metaData from response and stimulus fields into
    % subfields
    if isfield(packetCellArray{pp}.stimulus,'metaData')
        theConcatPacket.stimulus.metaData.(['packet' strtrim(num2str(pp))]) = ...
            packetCellArray{pp}.stimulus.metaData;
    end
    if isfield(packetCellArray{pp}.response,'metaData')
        theConcatPacket.response.metaData.(['packet' strtrim(num2str(pp))]) = ...
            packetCellArray{pp}.response.metaData;
    end
    if isfield(packetCellArray{pp},'metaData')
        theConcatPacket.metaData.(['packet' strtrim(num2str(pp))]) = ...
            packetCellArray{pp}.metaData;
    end
end % loop over packets

% Initialize the stimulus.values field, filled with the stimValueExtender
lengthStim=length(theConcatPacket.stimulus.timebase);
theConcatPacket.stimulus.values=zeros(totalStimValuesRows,lengthStim) + ...
    p.Results.stimValueExtender;

% Loop over packets and build the stimulus.values
stimRowCounter=1;
stimColCounter=1;
for pp=1:nPackets
    valuesRows=size(packetCellArray{pp}.stimulus.values,1);
    valuesCols=size(packetCellArray{pp}.stimulus.values,2);
    theConcatPacket.stimulus.values( ...
        stimRowCounter:stimRowCounter+valuesRows-1, ...
        stimColCounter:stimColCounter+valuesCols-1) = packetCellArray{pp}.stimulus.values;
    stimRowCounter=stimRowCounter+valuesRows;
    stimColCounter=stimColCounter+valuesCols;
end % loop over packets

% If a kernel field is defined, place this in theConcatPacket
if isfield(packetCellArray{1},'kernel')
    theConcatPacket.kernel=packetCellArray{1}.kernel;
end

end % function