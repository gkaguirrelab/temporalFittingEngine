function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the block temporal response model.
%
% Operates by calling forwardModelBTRM and then the applyKernel tfe method.
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
p.addParameter('addNoise',false,@islogical);

p.parse(params,stimulusStruct,kernelStruct,varargin{:});
params = p.Results.params;

%% Compute the forward model
modelResponseStruct = forwardModelBTRM(obj,params,stimulusStruct);

% report an iteration has completed
switch (obj.verbosity)
    case 'high'
        fprintf('.');
end

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add 1/f^2 noise
if p.Results.addNoise
    if ~isfield(params, 'noiseType')
        params.noiseType='white';
        
    end
switch params.noiseType
    case 'white'
        modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(modelResponseStruct.values));
    case 'red' % aka 'brown'ian noise or auto-regressive noise
        noise.timebase=modelResponseStruct.timebase;
        noise.values=normrnd(0,1,size(modelResponseStruct.values));
        noiseFilter.timebase=modelResponseStruct.timebase;
        noiseFilter.values=exp(-1/params.noiseTimeConstant*noiseFilter.timebase);
        noiseFilter=normalizeKernelArea(noiseFilter);
        noise=obj.applyKernel(noise,noiseFilter);
        noise.values=noise.values/std(noise.values)*params.noiseSd;
        modelResponseStruct.values = modelResponseStruct.values + noise.values;
    case 'none'
    otherwise
        error('That is not a noise type that I know');
end
end

end
