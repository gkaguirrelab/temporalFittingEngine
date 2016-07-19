function tmriDefaultParams(obj,varargin)

% Quadratic parameters
%
% We only store 5 because we handle the amplitude of the response
% in the amplitude parameter below.  The first axis of the ellipse
% has an implicit value of 1.
params.Qvec = [2 0.5 1 0 0]';

% Let's have a Naka-Rushton sigmoidal contrast response function
params.crfAmp = 1;
params.crfSemi = 1;
params.crfExponent = 2;

% Exponential falloff
params.expFalloff = 0.3;

% Pop structure into object
obj.params = params;

end