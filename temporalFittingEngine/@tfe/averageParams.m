function [meanParams,semParams] = averageParams(obj,paramsCellArray,varargin)
% plot(obj,timebase,response,varargin)
% Average tfe model parameters from multiple fits.
%
% Syntax:
%     [meanParams,semParams] = averageParams(obj,paramsCellArray);
%
% Description:
%
% Inputs:
%    paramsCellArray - Cell array where each entry is model parameters in
%                      the native format of the object.
%
% Outputs:
%    meanParams       - Mean params in the native format of the object.
%    semParams        - SEM of the parameters, also in the native format.
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

% Loop over all params to make a matrix of all params
for ii = 1:length(paramsCellArray)
    allParams(:,ii) = obj.paramsToVec(paramsCellArray{ii});
end

% Take the mean across params
meanVec = mean(allParams,2);

% Get the SEM
semVec = std(allParams,0,2)./sqrt(length(meanVec));

% Return params
meanParams = obj.vecToParams(meanVec);
semParams  = obj.vecToParams(semVec);

end