function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the block temporal response model.
%
% Operates by calling forwardModelLEAK and then the applyKernel tfe method.
%
% Optional key/value pairs
%   'AddNoise' - true/false (default false).  Add noise to computed
%     response? Useful for simulations.

%% Parse input.
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('params',@isstruct);
p.addRequired('stimulusStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) || isstruct(x)));
p.addParameter('addNoise','none',@ischar);
p.parse(params,stimulusStruct,kernelStruct,varargin{:});
params = p.Results.params;

%% Compute the forward model
modelResponseStruct = forwardModelLEAK(obj,params,stimulusStruct);

% report an iteration has completed
switch (obj.verbosity)
    case 'high'
        fprintf('.');
end

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add 1/f^2 noise
switch p.Results.addNoise
    case 'white'
        modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(modelResponseStruct.values));
    case 'red'
        noise.timebase=modelResponseStruct.timebase;
        noise.values=normrnd(0,params.noiseSd,size(modelResponseStruct.values));
        noiseFilter.timebase=modelResponseStruct.timebase;
        noiseFilter.values=exp(-1*params.noiseTau*noiseFilter.timebase);
        noiseFilter=normalizeKernelAmplitude(noiseFilter);
        noise=obj.applyKernel(noise,noiseFilter);
        modelResponseStruct.values = modelResponseStruct.values + noise.values;
    case 'none'
    otherwise
        error('That is not a noise addition method I know');
end


end
