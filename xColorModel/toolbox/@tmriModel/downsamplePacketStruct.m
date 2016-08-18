function downsampledStruct = downsamplePacketStruct(obj,S,downsampleFactor)

% function downsampledStruct = downsamplePacketStruct(S,downsampleFactor)
%
% downsamples timebase and values within a struct by a given factor
%
% inputs
%
% S               : struct with timebase and values fields
% downsampleFactor: factor by which to downsample

p = inputParser;
p.addRequired('S',@isstruct);
p.addRequired('downsampleFactor',@isnumeric);
p.parse(S,downsampleFactor);

% duplicate the input structure, but empty the timebase and values fields
downsampledStruct = S;
downsampledStruct.timebase = [];
downsampledStruct.values = [];

% loop over each row, and downsample
for i = 1:size(S.values,1)
    downsampledStruct.timebase(i,:) = downsample(S.timebase(i,:),downsampleFactor);
    downsampledStruct.values(i,:) = downsample(S.values(i,:),downsampleFactor);
end

end