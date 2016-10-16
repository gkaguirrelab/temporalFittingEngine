function [ theInstancePacket ] = splitOffAnInstancePacket( thePacket, params )
% function [ theInstancePacket ] = splitOffAnInstancePacket( thePacket, params )
%
% When passed a packet that contains multiple stimulus instances, the
% routine will split off a single instance (in both the stimulus and
% response fields) and return a packet that contains just this instance.
%
% The passed params argument defines the behavior of this split:
%   params.instanceIndex -- which instance to split off
%   params.splitOnsetMsecs -- the point to start the split
%   params.splitDurationMsecs -- the point to end the split
%   params.normFlag -- integer index flag, with the following behaviors:
%      1 (default) -- no normalization
%      2 -- recenter, original units
%      3 -- recenter, % change units
%   params.normalizationWindowMsecs -- Duration in msecs for the
%      normalization operations described in the prior two parameters. If
%      undefined, the normalizationWindow will be the entire response.
%


% check the passed params
if ~isfield(params,'instanceIndex')
    error('Need to indicate which instance to split off');
end
if ~isfield(params,'splitOnsetMsecs')
    % Find the stimulus onsets so that we can align the data to it. We
    % do that by finding a [0 1] edge from a difference operator.
    tmp = diff(thePacket.stimulus.values(params.instanceIndex,:));
    tmp(tmp < 0) = 0;
    tmp(tmp > 0) = 1;
    
    % Check if the very first value is 1, in which case the stim onset is
    % at the initial value
    if tmp(1)==1
        stimOnset = 1;
    else
        stimOnset = strfind(tmp, [0 1]);
    end
    
    if ~length(stimOnset)==1
        error('Cannot find a unique stimulus onset for this instance')
    end
    params.splitOnsetMsecs=stimOnset;
end
if ~isfield(params,'splitDurationMsecs')
    error('Need to indicate the extraction duration');
end
if ~isfield(params,'normFlag')
    params.normFlag=1;
end

% Map the input packet to the output packet
theInstancePacket=thePacket;

% Calculate the extraction indices for the stimulus vectors. We will first
% need the deltaT in the timebase of the stimulus
% derive the deltaT from the stimulusTimebase
check = diff(thePacket.stimulus.timebase);
deltaT = check(1);

% Find the index of the start of the extraction, rounding to the
% nearest index if necesary
[smallestDifference, startIndex] = min(abs(thePacket.stimulus.timebase - params.splitOnsetMsecs));

% Find the number of indices to extract, round to closest one.
nExtractionIndices=round(params.splitDurationMsecs/deltaT)-deltaT;

% Identify the stimulus indices to extract, making sure we are not
% exceeding the length of the stimulus vector
if (startIndex+nExtractionIndices) <= length(thePacket.stimulus.timebase)
    idxToExtract = startIndex:(startIndex+nExtractionIndices);
else
    idxToExtract = startIndex:length(thePacket.stimulus.timebase);
end

% Put the indices to extract back into the InstancePacket to be returned
theInstancePacket.stimulus.metaData.idxToExtract=idxToExtract;

% set theInstancePacket to have the extracted stimulus components
theInstancePacket.stimulus.timebase = thePacket.stimulus.timebase(idxToExtract);
theInstancePacket.stimulus.timebase = ...
    theInstancePacket.stimulus.timebase-thePacket.stimulus.timebase(idxToExtract(1));
theInstancePacket.stimulus.values = thePacket.stimulus.values(params.instanceIndex,idxToExtract);

% Calculate the extraction indices for the response vectors. We will first
% need the deltaT in the timebase of the stimulus
% derive the deltaT from the stimulusTimebase
check = diff(thePacket.response.timebase);
deltaT = check(1);

% Find the index of the start of the extraction, rounding to the
% nearest index if necesary
[smallestDifference, startIndex] = min(abs(thePacket.response.timebase - params.splitOnsetMsecs));

% Find the number of indices to extract, round to closest one.
nExtractionIndices=round(params.splitDurationMsecs/deltaT)-deltaT;

% Identify the stimulus indices to extract, making sure we are not
% exceeding the length of the stimulus vector
if (startIndex+nExtractionIndices) <= length(thePacket.response.timebase)
    idxToExtract = startIndex:(startIndex+nExtractionIndices);
else
    idxToExtract = startIndex:length(thePacket.response.timebase);
end

% set theInstancePacket to have the extracted response components
theInstancePacket.response.timebase=thePacket.response.timebase(idxToExtract);
theInstancePacket.response.timebase = ...
    theInstancePacket.response.timebase-thePacket.response.timebase(idxToExtract(1));
theInstancePacket.response.values=thePacket.response.values(idxToExtract);

%% Normalize the response time-series if requested

%   params.normFlag -- integer index flag, with the following behaviors:
%      1 (default) -- no normalization
%      2 -- recenter, original units
%      3 -- recenter, % change units
%   params.normalizationWindowMsecs -- Duration in msecs for the
%      normalization operations described in the prior two parameters. If
%      undefined, the normalizationWindow will be the entire response.

if ~isfield(params,'normalizationWindowMsecs')
    params.normalizationWindowMsecs=max(theInstancePacket.response.timebase);
end

% Set up the indices for normalization window
check = diff(theInstancePacket.response.timebase);
deltaT = check(1);
nNormIndices=round(params.normalizationWindowMsecs/deltaT)-deltaT;

% Shift the start of the norm window to the first, non-nan value in the
% response.values vector
nonNaNIndices=find(~isnan(theInstancePacket.response.values));
normStartIndex=nonNaNIndices(1);

% Get the normValue
normValue = ...
    nanmean(theInstancePacket.response.values(normStartIndex:normStartIndex+nNormIndices-1));

% Perform the normalization
switch params.normFlag
    case 1 % make no changes
        normValue=normValue;
    case 2 % recenter, original units
        theInstancePacket.response.values = ...
            theInstancePacket.response.values - normValue;
    case 3 % recenter, % change units
        theInstancePacket.response.values = ...
            (theInstancePacket.response.values - normValue)/normValue;
    otherwise
        error('Unknown normalization flag');
end


% place the parameters used for the extraction into the metadata field
theInstancePacket.metaData.splitOffAnInstance=params;

end

