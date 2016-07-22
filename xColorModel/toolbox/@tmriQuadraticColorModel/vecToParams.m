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

params.Qvec(1:5) = x(1:5)';
params.crfAmp = x(6);
params.crfExponent = x(7);
params.crfSemi = x(8);
params.expFalloff = x(9);

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    params.noiseSd = x(10);
end

end