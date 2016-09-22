function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the quadratic color model.
%
% Optional key/value pairs
%   'AddNoise' - true/false (default false).  Add noise to computed
%     response?  Useful for simulations.
%  'HRF' - structure (default empty).  Structure describing HRF to be used
%     to go from neural to BOLD response. If empty, no convolution is done

%% Parse input.
% At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('timebase',@isnumeric);
p.addRequired('stimulus',@isnumeric);
p.addParameter('AddNoise',false,@islogical);
p.addParameter('HRF',[]);
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
neuralResponse = path([params.crfAmp,params.crfSemi,params.crfExponent],theLengths);

%% Optionally, apply HRF
response = obj.applyKernal(timebase,neuralResponse,p.Results.HRF);

%% Optional add of noise
if (p.Results.AddNoise)
    response = response + normrnd(0,params.noiseSd,size(response));
end

end