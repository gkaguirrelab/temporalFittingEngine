% tmriQuadraticColorModelTest
%
% Test function for the quadratic color model.
%
% 6/26/16  dhb  Wrote it.

%% Clear and close
clear all; close all;

%% Construct the model object
tmri = tmriQuadraticColorModel;

%% Set parameters
%
% Six parameters define a unit quadratic form in three dimensions, I think.
params.Q = [1 0 0 1 0 0]';

% Let's have a Naka-Rushton sigmoidal contrast response function
params.crfAmp = 1;
params.crfExponent = 2;
params.crfSemi = 1;

% And an exponential falloff of response
params.expFalloff = 0.3;

% Tuck the parameter structure into the object
tmri.params = params;

%% Set the timebase we want to compute on
deltaT = 1;
totalTime = 1000;
timebase = 0:deltaT:totalTime;
tmri.timebase = timebase;

%% Specify the stimulus. 
%
% We'll specify this as a 3 by size(timebase,2) matrix,
% where each column is the signed L,M,S contrast of the stimulus
% at the specified time.
nTimeSamples = size(timebase,2);
stimulus = rand(3,nTimeSamples);
tmri.stimulus = stimulus;

%% Test that we can get a vector of paramters and put it back
x0 = tmri.paramstovec;
x1 = x0;
x1(2) = 2;
x1(7) = 3;
tmri.vectoparams(x1);
x2 = tmri.paramstovec;
if (any(x1 ~= x2))
    error('Parameter vectorizing and back not working right');
end
