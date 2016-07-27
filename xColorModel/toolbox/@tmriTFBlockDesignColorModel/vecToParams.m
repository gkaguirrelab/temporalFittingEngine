function params = vecToParams(obj,x,varargin)
% params = vecToParams(obj,x,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam'
%     true or false (default)

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('x',@isnumeric);
p.addParameter('UseNoiseParam',false,@islogical);
p.parse(x,varargin{:});
x = p.Results.x;

% Get base values of non vectorized parameters
params = obj.paramsBase;

% Push vector back into matrix in parameters structure, handling whether or
% not we had a noise parameter.
if (p.Results.UseNoiseParam)
    params.paramMainMatrix = reshape(x(1:end-1),params.matrixRows,params.matrixCols);
    params.noiseSd = x(end);
else
    params.paramMainMatrix = reshape(x,params.matrixRows,params.matrixCols);
end