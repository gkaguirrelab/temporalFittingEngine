function [paramsVec,vlbVec,vubVec] = tmriDefaultParams(obj,varargin)
% [params,vlb,vub] = tmriDefaultParams(obj,varargin)
%
% Set objects params to default values.
%
% Also return the default parameters in vector form as well as reasonable upper and lower
% bounds.

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

% Get vector form to return, as well as define bounds in vector form.
paramsVec = obj.paramsToVec;
vlbVec = [1e-3 1e-3 0 0 0 1e-3 1e-3 1e-2 1e-1]';
vubVec = [1e3 1e3 2*pi 2*pi 2*pi 1e3 1e3 1e2 1e1]';

end