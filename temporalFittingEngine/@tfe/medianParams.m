function [medianParams] = medianParams(obj,paramsCellArray,varargin)
% Take the median of tfe model parameters from multiple fits.
%
% Syntax:
%     [meanParams,semParams] = medianParams(obj,paramsCellArray);
%
% Description:
%
% Inputs:
%    paramsCellArray - Cell array where each entry is model parameters in
%                      the native format of the object.
%
% Outputs:
%    medianParams     - Median params in the native format of the object.
%
% Optional key/value pairs:
%    None.
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('paramsCellArray',@iscell);
p.parse(paramsCellArray,varargin{:});

% Check that paramsToVec method returns a column vector
checkVec = obj.paramsToVec(paramsCellArray{1});
if (~isvector(checkVec) | size(checkVec,2) > 1)
    error('paramsToVec method does not return a column vector');
end

% Loop over all params to make a matrix of all params
for ii = 1:length(paramsCellArray)
    allParams(:,ii) = obj.paramsToVec(paramsCellArray{ii});
end

% Take the mean across params
meanVec = nanmedian(allParams,2);

% Return params
medianParams = obj.vecToParams(meanVec);

end