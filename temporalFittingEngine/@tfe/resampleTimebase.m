function resampledStruct= resampleTimebase(obj,inputStruct,newTimebase,varargin)
% function resampledStruct= resampleTimebase(obj,inputStruct,newTimebase,varargin)
%
% Resamples timebase and values within an inputStruct to the passed
% new timebase.  Preserves the other fields of the struct.
%
% Required Inputs:
%   inputStruct: The input, must have timebase and values fields.
% 	  Typically this will be an elemement of a packet.
%   newTimebase: New timebase
%
% Time is always in milliseceonds.
%
% Optional key/value pairs:
%   'method' - string (default 'interp1_linear').  How to resample.
%     'interp1_linear' - Use Matlab's interp1, linear method and
%     extrapolate with 0 outside of input domain.

%% Parse input
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('inputStruct',@isstruct);
p.addRequired('newTimebase',@isnumeric);
p.addParameter('errorType', 'rmse', @ischar);
p.addParameter('method','interp1_linear',@ischar);
p.parse(inputStruct,newTimebase,varargin{:});

%% Resample
%
% Set up
resampledStruct = inputStruct;
nRows = size(inputStruct.values,1);
resampledStruct.timebase = newTimebase;
resampledStruct.values = zeros(nRows,length(newTimebase));

% Loop over each row of values, and downsample
for ii = 1:nRows
    switch (p.Results.method)
        case 'interp1_linear'
            resampledStruct.values(ii,:) = interp1(inputStruct.timebase,inputStruct.values(ii,:),newTimebase,'linear',0);
        otherwise
            error('Unknown method specified');
    end
end


end