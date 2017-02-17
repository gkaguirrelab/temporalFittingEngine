function [fVal,modelResponseStruct] = fitError(obj,paramsVec,thePacket,varargin)
% [fVal,modelResponseStruct] = fitError(obj,paramsVec,thePacket,varargin)
%
% Compute the error measure for passed model parameters, for the tfe class.
%
% Inputs:
%   paramsVec - model parameters in their vector form.  Just one of these
%   thePacket - standard packet containing stimulus, timebase, and metadata
%
% Optional key/value pairs
%  'errorType' - string (default 'rmse') Type of error to compute.
%    'rmse' - Root mean squared error.
%    '1-r2' - 1-R2
%
% Outputs: 
%   fVal: mean value of fit error, mean taken over runs.
%   modelResponseStruct: predicted response, standard structure form

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('thePacket',@isstruct);
p.addParameter('errorType','rmse',@ischar);
p.parse(paramsVec,thePacket,varargin{:});

%% Convert parameters to vector form
params = obj.vecToParams(paramsVec);

%% Compute the fit based on the timebase of the stimulus
modelResponseStruct = obj.computeResponse(params,thePacket.stimulus,thePacket.kernel,varargin{:});

%% Downsample to computed response to timebase of the response
modelResponseStruct = obj.resampleTimebase(modelResponseStruct,thePacket.response.timebase,varargin{:});

%% Get fit error measurement
switch (p.Results.errorType)
    case 'rmse'
        fVal = sqrt(nanmean((modelResponseStruct.values-thePacket.response.values).^2));
    case '1-r2'
        residuals=thePacket.response.values - modelResponseStruct.values;
        residuals=residuals-nanmean(residuals);
        fVal = nansum(residuals.^2)/nansum((thePacket.response.values - mean(thePacket.response.values)).^2);
    otherwise
        error('Unknown error type passed');
end

%% Check for inappropriate values in fVal. This is to aid debugging
if isnan(fVal)
    error('the model has returned an undefined error measurement');
end

end