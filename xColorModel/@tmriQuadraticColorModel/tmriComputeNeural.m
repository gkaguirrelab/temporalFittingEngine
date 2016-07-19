function obj = tmriComputeNeural(obj,varargin)
% obj = tmriComputeNeural(obj,varargin)
% 
% Compute method for the quadratic model. 

%% Get the ellipsoid parameters in cannonical form
Q = obj.Q;

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
theLengths = diag(obj.stimulus'*Q*obj.stimulus);

params.crfAmp = 1;
params.crfSemi = 1;
params.crfExponent = 2;

%% Push the quadratic response through a Naka-Rushton non-linearity
theResponse = ComputeNakaRushton([obj.params.crfAmp,obj.params.crfSemi,obj.params.crfExponent],theLengths);

%% Sock away the prediction
obj.neuralResponse = theResponse;

end