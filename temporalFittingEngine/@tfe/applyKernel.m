function convModelResponse = applyKernel(obj,timebase,modelResponse,theKernel,varargin)
% convModelResponse = applyKernel(obj,timebase,modelResponse,theKernel,varargin)
% 
% Apply a convolution kernel to the modeled response. In a typical
% application, this will be a hemodynamic response function applied to
% a model of neural activity to produce a BOLD fMRI response.
%
% Both modelResponse and theKernel arguments are structures, containing
% timebase and value fields.  The timebases do not need to be the same.
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
p.addRequired('theKernel',@(x)(isempty(x) || isstruct(x)));
p.parse(timebase,modelResponse,theKernel,varargin{:});

%% Convolve
%
% But if empty matrix is passed for kernel, do nothing,
if (isempty(theKernel) || isempty(theKernel.values))
    convModelResponse = modelResponse;
else 
    % Align kernel with 0, if not already done
    zeroAlignedKernel = theKernel.values-theKernel.values(1);
    
    % The kernel and modelResponse need to be same length: pad with 0's
    convKernel = zeros(size(modelResponse));
    convKernel(1:length(zeroAlignedKernel)) = zeroAlignedKernel;
    
    % Convolve stimulus with HRF to get BOLD response from neural response.
    convModelResponsePreCut = conv(modelResponse,convKernel) ;
    
    % Cut off extra conv values. Need to double check that this is right
    convModelResponse = convModelResponsePreCut(1:length(modelResponse));  
end

    
    