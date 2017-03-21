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

%% Optional add noise
if p.Results.addNoise
    if ~isfield(params, 'noiseInverseFrequencyPower')
        params.noiseInverseFrequencyPower=0;
    end
    cn = dsp.ColoredNoise(params.noiseInverseFrequencyPower,length(modelResponseStruct.timebase),1);
    noise = cn()';
    noise = noise/std(noise)*params.noiseSd;
    modelResponseStruct.values = modelResponseStruct.values + noise;
end

end
