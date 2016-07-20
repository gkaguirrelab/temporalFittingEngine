function response = tmriComputeNeural(obj,varargin)
% response = tmriComputeNeural(obj,varargin)
% 
% Compute method for the quadratic model. 

%% Get the ellipsoid parameters in cannonical form
Q = obj.simulateQ;
params = obj.simulateParams;

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
theLengths = diag(obj.stimulus'*Q*obj.stimulus);

%% Push the quadratic response through a Naka-Rushton non-linearity
theResponse = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],theLengths);

end