function [outputStruct,kernelStruct] = applyKernel(obj,inputStruct,kernelStruct,varargin)
% [outputStruct,kernelStruct] = applyKernel(obj,modelResponseStruct,kernelStruct,varargin)
%
% Apply a convolution kernel to a modeled response. In a typical
% application, this will be a hemodynamic response function applied to
% a model of neural activity to produce a BOLD fMRI response.
%
% The outputStruct's values field has the result of the convolution, and
% its timebase matches that of the input struct.
%
% The returend kernelStruct is the input, but if necessary its timebase has
% and values have been resampled to have the same delta time as the
% inputStruct's timebase.  The duration of the resampled kernel is at least
% as long as the original, and can be a little longer if the resampling
% requires an extension to produce an integer number of resampled times.
% This is returned mostly for debugging and checking purposes.
%
% Both modelResponseStruct and kernelStruct arguments are structures, containing
% timebase and values fields.  The timebases do not need to be the same, but
% each must be regularly sampled.
%
% Optional key/value pairs
%   'method' - string (default 'interp1_linear').  How to resample kernel timebase,
%      if needed.  This is passed onto method resampleTimebase.
%     'interp1_linear' - Use Matlab's interp1, linear method.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('inputStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) | isstruct(x)));
p.parse(inputStruct,kernelStruct,varargin{:});

%% Propagate all fields forward
outputStruct = inputStruct;

%% If empty matrix is passed for kernel, return
if (isempty(kernelStruct) || isempty(kernelStruct.values))
    return;
end

%% Get how many rows are in the the inputStruct value, and check
[nRows,nCols] = size(inputStruct.values);
if (size(inputStruct.timebase,2) ~= nCols)
    error('Badly formed response structure, length of timebase and values not the same.');
end
check = diff(inputStruct.timebase);
responseDeltaT = check(1);
if (any(abs(check - check(1)) > 1e-6))
    error('Response structure timebase is not regularly sampled');
end

%% Similar check on convolution kernel
if (length(kernelStruct.timebase) ~= length(kernelStruct.values))
    error('Badly formed kernel structure, length of timebase and values not the same.');
end
check = diff(kernelStruct.timebase);
kernelDeltaT = check(1);
if (any(abs(check - check(1)) > 1e-6))
    error('Kernel structure timebase is not regularly sampled');
end

%% Resample kernel to same delta time as response
if (responseDeltaT ~= kernelDeltaT)
    nSamples = ceil((kernelStruct.timebase(end)-kernelStruct.timebase(1))/responseDeltaT);
    newKernelTimebase = kernelStruct.timebase(1):responseDeltaT:(kernelStruct.timebase(1)+nSamples*responseDeltaT);
    kernelStruct = obj.resampleTimebase(kernelStruct,newKernelTimebase,varargin{:});
end

%% Loop over rows for the convolution
for ii=1:nRows
    % Convolve a row of inputStruct.values with the kernel.  The
    % convolutoin is a discrete approximation to an intergral, so we
    % explicitly include the factor of responseDeltaT.
    valuesRowConv = conv(inputStruct.values(ii,:),kernelStruct.values,'full')*responseDeltaT;
    
    % Cut off extra conv values
    outputStruct.values = valuesRowConv(1:length(inputStruct.timebase));  
end 

end


