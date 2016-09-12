function [fVal,responsePredicted] = fitError(obj,paramsVec,stimulus,responseToFit,HRF)
% [fVal,allFVals,responsePredicted] = fitError(obj,paramsVec,stimulus,responseToFit,HRF)
%
% Compute the error measure for passed model parameters, for the tmir
% class.
%
% Inputs:
%   paramsVec - model parameters in their vector form.  Just one of these
%   stimulus - struct containing stimulus, timebase, and metadata
%   responseToFit - struct containing response to fit
%   HRF - structure defining the HRF.  Can be empty, in which case no HRF
%         is applied.
%
% find the parameters that minimize the average fit error, taken over the runs.
% Doing it this way allows us to more easily use this routine for cross
% validation.
%
% Outputs: 
%   fVal: mean value of fit error, mean taken over runs.
%   responsePredicted: cell array of the predicted response, one for each run

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('stimulus',@isstruct);
p.addRequired('responseToFit',@isstruct);
p.addRequired('HRF',@isstruct);
p.parse(paramsVec,stimulus,responseToFit,HRF);

params = obj.vecToParams(paramsVec);

% compute the fit
responsePredictedPreDownsample = obj.computeResponse(params,stimulus.timebase,stimulus,'HRF',p.Results.HRF);
% downsample to resolution of data
responsePredicted = interp1(stimulus.timebase,responsePredictedPreDownsample,responseToFit.timebase);
% get error measurement
fVal = sqrt(mean((responsePredicted-responseToFit.values).^2));
display(num2str(fVal));

end