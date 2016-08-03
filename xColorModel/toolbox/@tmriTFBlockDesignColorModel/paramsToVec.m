function x = paramsToVec(obj,params,varargin)
% x = paramsToVec(obj,params,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam'
%     true or false (default) 

% Parse input
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('UseNoiseParam',false,@islogical);
p.parse(params,varargin{:});
params = p.Results.params;

% The parameters for search are in the paramMainMatrix field
x = reshape(params.paramMainMatrix,params.matrixRows*params.matrixCols,1);
obj.paramsBase = params;

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    x(end+1) = params.noiseSd;
end

end