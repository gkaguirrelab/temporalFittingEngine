function [concatParams] = concatenateParams(obj,paramsCellArray,varargin)
% Concatenate tfe model parameters from multiple fits.
%
% Syntax:
%     [concatPacket] = concatenateParams(obj,paramsCellArray,varargin)
%
% Description:
%
% Inputs:
%    paramsCellArray - Cell array where each entry is model parameters in
%                      the native format of the object.
%
% Outputs:
%    concatParams    - A packet containing the concatenated parameters.
%
% Optional key/value pairs:
%    averageBaseline  - averages the baseline condtition across the
%                       different paramters in the cell array. ASSUMES that 
%                       the baseline condtion is the last entry in each
%                       paramsMainMatrix subfield of the IAMP fit. 
%    makeBaselineZero - Make the baseline conditon 0. will be returned as
%                       the last entry of the params. 
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
p.addParameter('averageBaseline',false,@islogical);
p.addParameter('makeBaselineZero',true,@islogical);
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

% Handle baseline 
if p.Results.averageBaseline
    baseline = mean(allParams(end,:));
    allParams(end,:) = [];
elseif p.Results.makeBaselineZero
    baseline = 0;
    allParams(end,:) = [];
end

% Take the mean across params
concatVec = allParams(:);

if p.Results.averageBaseline || p.Results.makeBaselineZero
    concatVec = [concatVec;baseline];
end

% Return params
concatParams.paramNameCell   =  {'amplitude'};
concatParams.paramMainMatrix = concatVec';
concatParams.matrixRows      = size(concatVec,1);
concatParams.matrixCols      = size(concatVec,2);
concatParams.noiseSd         = 0;

end