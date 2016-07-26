% tmriTFBlockDesignColorModelDemo
%
% Demonstrate function for the quadratic color model.
%
% 6/26/16  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'..','BDCMToolbox'));
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'toolbox'));

%% Construct the model object
tmri = tmriTFBlockDesignColorModel;

%% Specify the stimulus
%
% Be sure to define defaultParamsInfo.nStimuli here, as we need it to create default
% parameters for the model.
defaultParamsInfo.nStimuli = 21;

%% Set parameters
%
% Six parameters define a quadratic form in three dimensions, but
% we normalize the first to 1 so we only need five numbers here.
params0 = tmri.defaultParams('DefaultParamsInfo',defaultParamsInfo);
fprintf('Default model parameters:\n');
%tmri.print(params0);

%% Test paramsToVec and vecToParams
params1 = params0;
params1.paramsMainMatrix = rand(size(params1.paramMainMatrix));
x1 = tmri.paramsToVec(params1);
params2 = tmri.vecToParams(x1);
if (any(params1.paramMainMatrix ~= params2.paramMainMatrix))
    error('vecToParams and paramsToVec do not invert');
end

%% Set the timebase we want to compute on
deltaT = 1;
totalTime = 1000;
timebase = 0:deltaT:totalTime;

%% Specify the stimulus. 
%
% We'll specify this as a 3 by size(timebase,2) matrix,
% where each column is the signed L,M,S contrast of the stimulus
% at the specified time.  And then we'll blur it so that we have
% a smoothish signal to look at.
nTimeSamples = size(timebase,2);
filter = fspecial('gaussian',[1 nTimeSamples],6);
stimulus= rand(3,nTimeSamples);
for i = 1:3
    % stimulus(i,:) = conv(stimulus(i,:),filter,'same');
    stimulus(i,:) = ifft(fft(stimulus(i,:)) .* fft(filter)); 
end

%% Test that we can obtain a neural response
params1.crfAmp = 2;
params1.crfSemi = 0.5;
params1.crfExponent = 3;
params1.noiseSd = 0.02;
fprintf('Simulated model parameters:\n');
tmri.print(params1);
responseToFit = tmri.computeResponse(params1,timebase,stimulus,'AddNoise',true);
tmri.plot(timebase,responseToFit);

%% Test the fitter
[paramsFit,fitResponse] = tmri.fitResponse(timebase,stimulus,responseToFit,'DefaultParamsInfo',defaultParamsInfo);
fprintf('Model parameter from fits:\n');
tmri.print(paramsFit);
tmri.plot(timebase,fitResponse,'Color',[0 1 0],'NewWindow',false);

