function [outputStruct,kernelStruct] = applyKernel(obj,inputStruct,kernelStruct,varargin)
% Apply a convolution kernel to a passed inputStruct
%
% Syntax:
%  [outputStruct,kernelStruct] = applyKernel(obj,modelResponseStruct,kernelStruct,varargin)
%
% Description:
%	Apply a convolution kernel to a modeled response. In a typical
%   application, this will be a hemodynamic response function applied to a
%   model of neural activity to produce a BOLD fMRI response.
%
%   The outputStruct's values field has the result of the convolution, and
%   its timebase matches that of the input struct.
%
%   The returned kernelStruct is the input, but if necessary its timebase
%   has and values have been resampled to have the same delta time as the
%   inputStruct's timebase.  The duration of the resampled kernel is at
%   least as long as the original, and can be a little longer if the
%   resampling requires an extension to produce an integer number of
%   resampled times. This is returned mostly for debugging and checking
%   purposes.
%
%   Both modelResponseStruct and kernelStruct arguments are structures,
%   containing timebase and values fields.  The timebases do not need to be
%   the same, but each must be regularly sampled.
%
% Inputs
%   inputStruct           - F
%   kernelStruct          - F
%
% Optional key/value pairs
%  'method'               - String (default 'interp1_linear').  How to
%                           resample kernel timebase, if needed. This is
%                           passed onto method resampleTimebase.
%  'interp1_linear'       - Use Matlab's interp1, linear method.
%
% Outputs
%   outputStruct          - F
%   kernelStruct          - F
%
% Examples:
%{
    % Standard convolution. Test result against cached value.
    % We create a delta function response and a double-gamma HRF
    responseStruct.timebase = 0:1:25999;
	responseStruct.values = zeros(1,length(responseStruct.timebase));
	responseStruct.values(2000) = 1;
    kernelStruct.timebase=linspace(0,15999,16000);
    hrf = gampdf(kernelStruct.timebase/1000, 6, 1) -  gampdf(kernelStruct.timebase/1000, 12, 1)/10;
    kernelStruct.values=hrf;
    [ kernelStruct ] = normalizeKernelArea( kernelStruct );
    % Instantiate the tfe and perform the convolution
    temporalFit = tfeIAMP('verbosity','none');
    convResponseStruct = temporalFit.applyKernel(responseStruct,kernelStruct);
    % Compare the output to a cached hash of the output
    cachedHash = '41d0741a71e625ecc91c67f227017425';
    computedHash = DataHash(convResponseStruct);
    assert(strcmp(computedHash, cachedHash));
%}


%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('inputStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) | isstruct(x)));
p.addRequired('durationMsecsOfNansToCensor',5000,@isscalar);
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

%% Find the time of the absolute maximum of the kernel
% This may be used if we are handling nan values in the input below
% DO THIS HERE


%% Loop over rows for the convolution
for ii=1:nRows

    % Grab this row of the input structure
    inputRow = inputStruct.values(ii,:);

    % Detect if nans are present in this row of values
    nanIdx = isnan(inputRow);
    if ~isempty(nanIdx)
        % There are nans present. Spline interpolate over the missing
        % values
        inputRow=spline(inputStruct.timebase(~nanIdx), inputRow(~nanIdx), inputStruct.timebase);
    end

    % Convolve with the kernel.  The convolution is a discrete
    % approximation to an intergral, so we explicitly include the factor of
    % responseDeltaT.
    valuesRowConv = conv(inputRow,kernelStruct.values, 'full')*responseDeltaT;

    % Cut off extra conv values
    outputStruct.values(ii,:) = valuesRowConv(1:length(inputStruct.timebase));
    
    % If nans were present in the input vector, find those stretches of
    % nans that are longer than the passed threshold and set the output
    % vector to have nans in the corresponding, temporally shifted
    % location.
    if ~isempty(nanIdx)
        % Find the stretches
        % nan the outputStruct.values at the corresponding delayed
        % locations
    end

    %     % new nan-safe method
    %     if strcmp(p.Results.convolveMethod, 'nanconv')
    %         outputStruct.values(ii,:) = nanconv_local(inputStruct.values(ii,:),kernelStruct.values, '1d')*responseDeltaT;
    %         % note that this function already cuts off the extra convolve
    %         % values, so we don't need to do it after
    %     end
    %
    %     % principled method: If the missle values are within the temporal
    %     % domain of the kernel (the length of time that the kernel spans), we
    %     % will interpolate the missing values and continue with the convolution
    %     % as normal. If, however, the length of the stretch of NaN values is
    %     % greater than the domain of the kernel, we will assign the
    %     % corresponding region of the output vector to be NaN as well.
    %
    %     if strcmp(p.Results.convolveMethod, 'principled')
    %         % identify any NaN values
    %         NaNIndices = find(isnan(inputStruct.values(ii,:)));
    %
    %         % determine the length of each NaN segment
    %
    %         % interpolate over NaN segments that are short enough in duration,
    %         % relative to kernel temporal domain
    %
    %         % lingering question: how to handle the convolution when we will
    %         % not be interpolating
    %         % 1st option: interpolate, then convolve, then replace these indices
    %         % with NaN
    %         % 2nd option: keep NaN, use nanconv, then replace these indices
    %         % with NaN.
    %         % These two options will produce different results for values along
    %         % the edges of these NaN segments
    %
    %     end
    
    
    
end



