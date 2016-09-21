function convModelResponse = applyKernel(obj,timebase,modelResponse,theKernel,varargin)
% boldResponse = applyKernel(obj,timebase,modelResponse,theKernel,varargin)
% 
% Apply a convolution kernel to the modeled response. In a typical
% application, this will be a hemodynamic response function applied to
% a model of neural activity to produce a BOLD fMRI response.
%
% We assume that the kernel structure contains a single field, kernel, on
% the same timebase spacing as the model response.
%
% Inputs:
%   timebase - times on which data/model predictions exist (in msecs)
%   modelResponse - the modeled response on timebase
%   theKernel - structure containing info about the convolution kernel
% 
% Optional key/value pairs

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('timebase',@isnumeric);
p.addRequired('modelResponse',@isnumeric);
p.addRequired('theKernel',@isstruct);
p.parse(timebase,modelResponse,theKernel,varargin{:});

%% If empty matrix is passed for kernel, do nothing,
if (isempty(theKernel.values))
    convModelResponse = modelResponse;
else %   convolve with kernel.
    % align kernel with 0, if not already done
    zeroAlignedKernel = theKernel.values-theKernel.values(1);
    % the kernel and modelResponse need to be same length: pad with 0's
    convKernel = zeros(size(modelResponse));
    convKernel(1:length(zeroAlignedKernel)) = zeroAlignedKernel;
    % Convolve stimulus with HRF to get BOLD response from neural response.
    convModelResponsePreCut = conv(modelResponse,convKernel) ;
    
    % Cut off extra conv values
    % Need to double check 
    convModelResponse = convModelResponsePreCut(1:length(modelResponse));  
end
