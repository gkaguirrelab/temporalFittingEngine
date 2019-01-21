function [concatParams,baselineShift] = concatenateParams(obj,paramsCellArray,varargin)
% Concatenate tfe model parameters from multiple fits.
%
% Syntax:
%     [concatPacket] = concatenateParams(obj,paramsCellArray,varargin)
%
% Description:
%     This concatenates together data from two IAMP fits, By default,
%     just concatenates, but you can control the behavior for the conventions
%     of some of our experiments by key/value pairs as described below.
%
%     A bit specifiç to our contrast/response experiment, which is
%     to say that someone could run an IAMP fit where the last parameter
%     didn't correspond to a baseline measurment.  You can override this
%     behavior with a key/value pair.
%
%     The noiseSd field of the returned parameters, which is only used in
%     forward simulation, is set to whatever it was in the first passed set
%     of parameters.
%
% Inputs:
%    paramsCellArray - Cell array where each entry is model parameters in
%                      the native format of the object.
%
% Outputs:
%    concatParams    - A packet containing the concatenated parameters.
%    baselineShift   - This contains the value we adjusted each sessions'
%                      baseline by.  Add this back in to predictions to
%                      match original time course data appropriatetly.  If
%                      no baseline handling was specified by a key/value
%                      pair, this is a vector of zeros.
%
% Optional key/value pairs:
%    'baselineMethod' - Baseline handling as below.  Default is 'none'.
%                         - 'none' - Don't do anything about baseline, just
%                         concatenate all parameters. - 'averageBaseline'
%                         - averages the baseline condtition across the
%                         different paramters in the cell array. ASSUMES
%                         that the baseline condtion is the last entry in
%                         each paramsMainMatrix subfield of the IAMP fit.
%                         Leaves other parameters alone and baselineShift
%                         vector is returned as zeros because when we do
%                         this we're viewing the variability of the
%                         baseline as noise.
%                         - 'makeBaselineZero' - Make the baseline conditon
%                         0. will be returned as the last entry of the
%                         params, and the other parameters are adjusted by
%                         subtracting off the original baseline parameter
%                         for each session.
%

% History:
%   01/14/19 mab 

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('paramsCellArray',@iscell);
p.addParameter('baselineMethod','none',@ischar);
p.parse(paramsCellArray,varargin{:});

% Check that paramsToVec method returns a column vector
checkVec = obj.paramsToVec(paramsCellArray{1});
if (~isvector(checkVec) | size(checkVec,2) > 1)
    error('paramsToVec method does not return a column vector');
end

% Check dimensions of passed cell array
if (size(paramsCellArray,1) > 1 & size(paramsCellArray,2) > 1)
    error('Passed cell array of parameters is in matrix form');
end

% Loop over all params to make a matrix of all params
for ii = 1:length(paramsCellArray)
    allParams(:,ii) = obj.paramsToVec(paramsCellArray{ii});
    baselineValues(1,ii) = allParams(end,ii);
end

% Handle baseline 
baselineShift = zeros(length(paramsCellArray),1);
switch (p.Results.baselineMethod)
    case 'averageBaseline'
        baseline = mean(allParams(end,:));
        allParams(end,:) = [];
    case 'makeBaselineZero'
        baseline = 0;
        allParams(end,:) = [];
        allParams = bsxfun(@minus,allParams,baselineValues);
        baselineShift = baselineValues';
    case 'none'
    otherwise
        error('Unknown baseline method passed');
end

% Take the mean across params
concatVec = allParams(:);

% Actually adjust baseline if we need to.
switch (p.Results.baselineMethod)
    case {'averageBaseline', 'makeBaselineZero'}
        concatVec = [concatVec ; baseline];  
end

% Return params
concatParams.paramNameCell   =  {'amplitude'};
concatParams.paramMainMatrix = concatVec';
concatParams.matrixRows      = size(concatVec,1);
concatParams.matrixCols      = size(concatVec,2);
concatParams.noiseSd         = paramsCellArray{1}.noiseSd;

end