function [fVal,responsePredicted] = fitError(obj,paramsVec,stimulus,responseToFit,kernel)
% [fVal,allFVals,responsePredicted] = fitError(obj,paramsVec,stimulus,responseToFit,HRF)
%
% Compute the error measure for passed model parameters, for the tfe class.
%
% Inputs:
%   paramsVec - model parameters in their vector form.  Just one of these
%   stimulus - struct containing stimulus, timebase, and metadata
%   responseToFit - struct containing response to fit
%   kernel - structure defining the kernel.  Can be empty, in which case no kernel
%         is applied.
%
% Optional key/value pairs
%  'errorType' - string (default 'rmse') Type of error to compute.
%    'rmse' - Root mean squared error.
%
% Outputs: 
%   fVal: mean value of fit error, mean taken over runs.
%   responsePredicted: cell array of the predicted response, one for each run

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('stimulus',@isstruct);
p.addRequired('responseToFit',@isstruct);
p.addRequired('kernel',@isstruct);
p.parse(paramsVec,stimulus,responseToFit,kernel);

%% Convert parameters to vector form
params = obj.vecToParams(paramsVec);

%% Compute the fit based on the timebase of the stimulus
responsePredictedPreDownsample = obj.computeResponse(params,stimulus.timebase,stimulus,'kernel',p.Results.kernel);

%% Downsample to computed response to timebase of the response
responsePredicted = interp1(stimulus.timebase,responsePredictedPreDownsample,responseToFit.timebase);

%% Get fit error measurement
fVal = sqrt(mean((responsePredicted-responseToFit.values).^2));

end