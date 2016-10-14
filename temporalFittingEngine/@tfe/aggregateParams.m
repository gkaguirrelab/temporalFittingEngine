function [centralParams] = aggregateParams(obj,paramsCellArray,varargin)
% function [paramsFitCentralTendency] = combineParamsFit(obj,paramsFitCellArray,varargin)
%
% Computes the central tendency of a set of passed paramsFit variabled
% included in the cell array.
%
% Optional key/value pairs
%   'combineMethod' - string (default 'mean'). How to combine the params.
%       mean
%       median
%       euclidean
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('paramsCellArray',@iscell);
p.addParameter('aggregateMethod','mean',@ischar);
p.parse(paramsCellArray,varargin{:});

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
    case 'median'
        centralParamsMainMatrix = median(allParamMainMatrix,1);
    case 'euclidean'
        centralParamsMainMatrix = sqrt(sum(allParamMainMatrix.^2,1));
    otherwise
        error('Unknown central tendency method passed');
end

centralParams = paramsCellArray{1};
centralParams.paramMainMatrix = centralParamsMainMatrix;

end


