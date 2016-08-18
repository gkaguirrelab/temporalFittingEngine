function boldResponse = applyHRF(obj,timebase,neuralResponse,theHRF,varargin)
% boldResponse = applyHRF(obj,timebase,neuralResponse,theHRF,varargin)
% 
% Apply the HRF (contained in the theHRF structure) to the neural response to
% obtain the BOLD response.  For right now, we'll assume that the HRF
% structure contains a single field, HRF, on the same timebase spacing as
% the neural response.
%
% Inputs:
%   timebase - times on which data/model predictions exist (in seconds)
%   neuralResponse - neural response on timebase
%   theHRF - structure containing info about the HRF
% 
% Optional key/value pairs

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('timebase',@isnumeric);
p.addRequired('neuralResponse',@isnumeric);
p.addRequired('HRF',@isstruct);
p.parse(timebase,neuralResponse,theHRF,varargin{:});

%% If empty matrix is passed for HRF, do nothing, otherwise convolve with HRF.
if (isempty(theHRF))
    boldResponse = neuralResponse;
else 
    % align HRF with 0, if not already done
    zeroAlignedHRF = theHRF.values-theHRF.values(1);
    % HRF and neuralResponse need to be same length: pad with 0's
    HRFconvKernel = zeros(size(neuralResponse));
    HRFconvKernel(1:length(zeroAlignedHRF)) = zeroAlignedHRF;
    % Convolve stimulus with HRF to get BOLD response from neural response.
    boldResponsePreCut = conv(neuralResponse,HRFconvKernel) ;
    
    % Cut off extra conv values
    %
    % Need to double check 
    boldResponse = boldResponsePreCut(1:length(neuralResponse));  
end
