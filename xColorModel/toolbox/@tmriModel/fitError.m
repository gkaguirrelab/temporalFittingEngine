function [fVal,allFVals,responsePredicted] = fitError(obj,paramsVec,timebase,stimulus,responseToFit,HRF)
% [fVal,allFVals,responsePredicted] = fitError(obj,paramsVec,timebase,stimulus,responseToFit,HRF)
%
% Compute the error measure for passed model parameters, for the tmir
% class.
%
% Inputs:
%   paramsVec - model parameters in their vector form.  Just one of these
%   timebase - cell array of times on which data/model predictions exist
%   stimulus - cell array of stimulus descriptions
%   responseToFit - cell array of responses to fit
%   HRF - structure defining the HRF.  Can be empty, in which case no HRF
%         is applied.
%
% Each entry of the cell arrays is for one run, and the fit is to find the
% parameters that minimize the average fit error, taken over the runs.
% Doing it this way allows us to more easily use this routine for cross
% validation.
%
% Outputs: 
%   fVal: mean value of fit error, mean taken over runs.
%   allFVals: vector of fit errors for the individual runs.
%   predictedResponse: cell array of the predicted response, one for each run

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('timebase',@iscell);
p.addRequired('stimulus',@iscell);
p.addRequired('responseToFit',@iscell);
p.addRequired('HRF');
p.parse(paramsVec,timebase,stimulus,responseToFit,HRF);

params = obj.vecToParams(paramsVec);
for ii = 1:length(timebase)
    responsePredicted{ii} = obj.computeResponse(params,timebase{ii},stimulus{ii},'HRF',p.Results.HRF);
    allFVals(ii) = sqrt(mean((responsePredicted{ii}-responseToFit{ii}).^2));
end
fVal = mean(allFVals);

end