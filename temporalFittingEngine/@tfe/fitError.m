function [fVal,modelResponseStruct] = fitError(obj,paramsVec,thePacket,varargin)
% [fVal,modelResponseStruct] = fitError(obj,paramsVec,thePacket,varargin)
%
% Compute the error measure for passed model parameters, for the tfe class.
%
% Inputs:
%   paramsVec           - Model parameters in their vector form.  Just one of these
%   thePacket           - Standard packet containing stimulus, timebase, and metadata
%
% Outputs: 
%   fVal                - Fit error
%   modelResponseStruct - Predicted response, standard structure form
%
% Optional key/value pairs
%  'errorType'         - String (default 'rmse') Type of error to compute.
%                         'rmse' - Root mean squared error.
%                         '1-r2' - 1-r2
%  'errorWeightVector' - Vector of weights to use on error for each
%                        response value. Only valid if error type is 'rmse'.
%                        Default empty.
%  'fitErrorScalar'    - Computed fit error is multiplied by this before
%                        return.  Sometimes getting the objective function
%                        onto the right scale makes all the difference in
%                        fitting. Default 1.
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('thePacket',@isstruct);
p.addParameter('errorType','rmse',@ischar);
p.addParameter('errorWeightVector',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fitErrorScalar',1,@isnumeric);
p.parse(paramsVec,thePacket,varargin{:});

%% If passed, perform some errorWeightVector checks
if ~isempty(p.Results.errorWeightVector)
    if length(p.Results.errorWeightVector)~=length(thePacket.response.values)
        error('The errorWeightVector is not the same length as the response');
    end
    if ~strcmp(p.Results.errorType,'rmse')
        error('An errorWeightVector can only be used with an rmse errorType');
    end
end

%% Convert parameters to vector form
params = obj.vecToParams(paramsVec);

%% Compute the fit based on the timebase of the stimulus
modelResponseStruct = obj.computeResponse(params,thePacket.stimulus,thePacket.kernel,varargin{:});

%% Downsample to computed response to timebase of the response
modelResponseStruct = obj.resampleTimebase(modelResponseStruct,thePacket.response.timebase,varargin{:});

%% Get fit error measurement
residuals=thePacket.response.values - modelResponseStruct.values;
switch (p.Results.errorType)
    case 'rmse'
        errorVector=residuals.^2;
        if ~isempty(p.Results.errorWeightVector)
            errorVector=errorVector*p.Results.errorWeightVector;
        end
        fVal = sqrt(nanmean(errorVector));
    case '1-r2'
        fVal = nansum((residuals-nanmean(residuals)).^2)/nansum((thePacket.response.values - mean(thePacket.response.values)).^2);
    otherwise
        error('Unknown error type passed');
end

%% Check for inappropriate values in fVal. This is to aid debugging
if isnan(fVal)
    error('the model has returned an undefined error measurement');
end

%% Multiply by fitErrorScalar
fVal = p.Results.fitErrorScalar*fVal;

end