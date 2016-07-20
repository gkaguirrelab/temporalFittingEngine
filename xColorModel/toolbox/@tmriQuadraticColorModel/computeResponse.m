function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the quadratic model. 

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('timebase',@isnumeric);
p.addRequired('stimulus',@isnumeric);
p.addParameter('PrintType','parameters',@isstring);
p.parse(params,timebase,stimulus,varargin{:});
params = p.Results.params;
timebase = p.Results.timebase;
stimulus = p.Results.stimulus;

%% Get the ellipsoid parameters in cannonical form
[~,~,Q] = EllipsoidMatricesGenerate([1 params.Qvec]');

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
theLengths = diag(stimulus'*Q*stimulus);

%% Push the quadratic response through a Naka-Rushton non-linearity
response = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],theLengths);

end