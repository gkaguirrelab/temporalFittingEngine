function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct form, use paramsToVec on the structs to
% get vector form.

%% Default parameters
param.startTime = 100; % time (in msecs) that the initial transient pupil response begins

% parameters of the initial gamma convolution
param.gammaTau = 200;    % time constant of the Gamma function (msecs)

% parameters of the sustained response
param.sustainedAmp = 0.2; % scaling of the sustained component
param.sustainedTau = 500;  % time constant of the low-pass (exponential decay) component (msecs).

% parameters of the persistent
param.persistentAmp = 0.2; % amplitude scaling of the persistent response
param.persistentT50 = 250; % time to half-peak of the super-saturating function (msecs)
param.persistentAlpha = 1500;  % time constant of the decay of the super-saturating function (msecs).

%% Lower bounds
paramsLb.startTime = 100;
paramsLb.gammaTau = 100;
paramsLb.sustainedAmp = 0;
paramsLb.sustainedTau = 0;
paramsLb.persistentAmp = 0;
paramsLb.persistentT50 = 0;
paramsLb.persistentAlpha = 500;

%% Upper bounds
paramsUb.startTime = 400;
paramsUb.gammaTau = 300;
paramsUb.sustainedAmp = 1;
paramsUb.sustainedTau = 2000;
paramsUb.persistentAmp = 1;
paramsUb.persistentT50 = 2000;
paramsUb.persistentAlpha = 6000;

end