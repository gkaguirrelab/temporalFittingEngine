function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the quadratic color model.
%
% Optional key/value pairs
%   'AddNoise' - true/false (default false).  Add noise to computed
%     response?  Useful for simulations.

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('stimulusStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) || isstruct(x)));
p.addParameter('addNoise',false,@islogical);
p.parse(params,stimulusStruct,kernelStruct,varargin{:});
params = p.Results.params;

%% Get the ellipsoid parameters in cannonical form
[~,~,Q] = EllipsoidMatricesGenerate([1 params.Qvec]');

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
theLengths = diag(stimulusStruct.values'*Q*stimulusStruct.values);

%% Push the quadratic response through a Naka-Rushton non-linearity
neuralResponse = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],theLengths);

%% Make the neural response structure
modelResponseStruct.timebase = stimulusStruct.timebase;
modelResponseStruct.values = neuralResponse;

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add noise
if (p.Results.addNoise)
    modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(modelResponseStruct.values));
end

end
