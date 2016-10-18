function [centralParams, centralfVal] = aggregateParams(obj,paramsCellArray, fValsArray, varargin)
% function [centralParams] = aggregateParams(obj,paramsCellArray,varargin)
%
% Computes the central tendency of a set of passed paramsFit variables
% included in the cell array.
%
% Optional key/value pairs
%   'aggregateMethod' - string (default 'mean'). How to combine the params.
%       mean
%       median
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('paramsCellArray',@iscell);
p.addRequired('fValsArray',@isnumeric);
p.addParameter('aggregateMethod','mean',@ischar);
p.parse(paramsCellArray,varargin{:});


%% Check that the items in the param cell array have the same parameter
% names in the same order

% Setup
nParamSets = length(paramsCellArray);

% Loop over the elements of the cell Array and build a matrix of params
for ss=1:nParamSets
    allParamMainMatrix(ss,:,:)=paramsCellArray{ss}.paramMainMatrix;
end

% Take the central tendency according to the specified method 
switch (p.Results.aggregateMethod)
    case 'mean'
        centralParamsMainMatrix = mean(allParamMainMatrix,1);
        centralfVal = mean(fValsArray);
    case 'median'
        centralParamsMainMatrix = median(allParamMainMatrix,1);
        centralfVal = median(fValsArray);
    otherwise
        error('Unknown aggregation method requested');
end

centralParams = paramsCellArray{1};
centralParams.paramMainMatrix = centralParamsMainMatrix;

end


