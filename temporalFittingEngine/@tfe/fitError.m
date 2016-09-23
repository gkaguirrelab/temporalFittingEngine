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
%
% Outputs: 
%   fVal: mean value of fit error, mean taken over runs.
%   modelResponsStruct: predicted response, standard structure form

%% Parse vargin for options passed here
p = inputParser;
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
        fVal = sqrt(mean((modelResponseStruct.values-thePacket.response.values).^2));
    otherwise
        error('Unknown error type passed');
end

end