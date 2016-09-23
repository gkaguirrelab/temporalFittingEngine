function convModelResponseStruct = applyKernel(obj,modelResponseStruct,kernelStruct,varargin)
% convModelResponseStruct = applyKernel(obj,modelResponseStruct,kernelStruct,varargin)
%
% THIS NEEDS TO BE WRITTEN CAREFULLY AND HAVE A GOOD UNIT TEST
% MAKE SURE TO DEAL WITH INCONSISTENT TIMEBASES
% 
% Apply a convolution kernel to a modeled response. In a typical
% application, this will be a hemodynamic response function applied to
% a model of neural activity to produce a BOLD fMRI response.
%
% Both modelResponseStruct and kernelStruct arguments are structures, containing
% timebase and value fields.  The timebases do not need to be the same.
% 
% Optional key/value pairs

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('modelResponseStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isemtpy(x) | isstruct(x)));
p.parse(modelResponseStruct,kernelStruct,varargin{:});

%% Propagate all fields forward
convModelResponseStruct = modelResponseStruct;

%% Convolve
%
% But if empty matrix is passed for kernel, do nothing,
if (isempty(kernelStruct) || isempty(kernelStruct.values))
else  
    % Align kernel with 0, if not already done
    zeroAlignedKernel = kernelStruct.values-kernelStruct.values(1);
    
    % The kernel and modelResponseStruct need to be same length: pad with 0's
    convKernel = zeros(size(modelResponseStruct.values));
    convKernel(1:length(zeroAlignedKernel)) = zeroAlignedKernel;
    
    % Convolve stimulus with HRF to get BOLD response from neural response.
    convModelResponsePreCut = conv(modelResponseStruct.values,convKernel) ;
    
    % Cut off extra conv values. Need to double check that this is right
    convModelResponseStruct.values = convModelResponsePreCut(1:length(modelResponseStruct.values));  
end

    
    