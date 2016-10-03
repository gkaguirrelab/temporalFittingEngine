function resampledStruct= resampleTimebase(obj,inputStruct,newTimebase,varargin)
% function resampledStruct= resampleTimebase(obj,inputStruct,newTimebase,varargin)
%
% Resamples timebase and values within an inputStruct to the passed
% new timebase.  Preserves the other fields of the struct.
%
% Required Inputs:
%   inputStruct: The input, must have timebase and values fields.
% 	  Typically this will be a structure taken from a packet.
%   newTimebase: New timebase
%
% Optional key/value pairs:
%   'method' - string (default 'linear'). How to resample.
%     'interp1_linear' - Use Matlab's resample linear method and
%          extrapolate with 0 outside of input domain.
%     'resample - Use Matlab's resample function

%% Parse input
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('inputStruct',@isstruct);
p.addRequired('newTimebase',@isnumeric);
p.addParameter('errorType', 'rmse', @ischar);
p.addParameter('resampleMethod','interp1_linear',@ischar);
p.parse(inputStruct,newTimebase,varargin{:});

% Transer all fields from the inputStruct to the resampledStruct
resampledStruct = inputStruct;

% How rows are there in the values field to resample?
nRows = size(inputStruct.values,1);

% Set up the resampledStruct timebase and values
resampledStruct.timebase = newTimebase;
resampledStruct.values = zeros(nRows,length(newTimebase));

% Loop over each row of values, and downsample
for ii = 1:nRows
    switch (p.Results.resampleMethod)
        case 'interp1_linear'
            resampledStruct.values(ii,:) = interp1(inputStruct.timebase,inputStruct.values(ii,:),newTimebase,'linear',0);
        case 'resample'
            % create a Matlab timeseries object out of elements of the
            % inputStruct
            inputRowTimeSeriesObject = ...
                timeseries(inputStruct.values(ii,:),inputStruct.timebase);
            resampledRowTimeSeriesObject = ...
                resample(inputRowTimeSeriesObject,newTimebase);
            resampledStruct.values(ii,:)=resampledRowTimeSeriesObject.data;
        otherwise
            error('Unknown method specified');
    end
end


end