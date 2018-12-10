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

    % A convolution containing varying amounts of NaNs
    responseStruct.timebase = 0:1:25999;
    % create a sine wave responseStruct
    sinFunc = @sin;
    responseStruct.values = sinFunc(responseStruct.timebase/1000);
    % add some NaN values
    responseStruct.values(10:20) = NaN;
    responseStruct.values(50:620) = NaN;
    % note the chunk that follows will be longer than the specified
    % 'durationMsecsOfNansToCensor' and will ultimately be censored
    responseStruct.values(10000:16000) = NaN; 
    % create standard HRF
    kernelStruct.timebase=linspace(0,15999,16000);
    hrf = gampdf(kernelStruct.timebase/1000, 6, 1) -  gampdf(kernelStruct.timebase/1000, 12, 1)/10;
    kernelStruct.values=hrf;
    [ kernelStruct ] = normalizeKernelArea( kernelStruct );
    % Instantiate the tfe and perform the convolution
    temporalFit = tfeIAMP('verbosity','none');
    convResponseStruct = temporalFit.applyKernel(responseStruct,kernelStruct, 'durationMsecsOfNansToCensor', 5000);
    figure; hold on;
    plot(responseStruct.timebase, responseStruct.values);
    plot(convResponseStruct.timebase, convResponseStruct.values);
    legend('Original Response', 'Convolved Response')
    xlabel('Time')
    ylabel('Response')
    % on this plot, notice that the first couple of NaN runs have been interpolated
    and their corresponding values are present in the final convolution
    product. The longer NaN run, however, is ultimately censored in the
    final output, after introducing the proper timeshift.
%}


%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('inputStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) | isstruct(x)));
p.addParameter('durationMsecsOfNansToCensor',5000,@isscalar);
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
[lagToPeakValue, lagToPeakIndex] = max(abs(kernelStruct.values));


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
    nanRuns = [];
    if ~isempty(find(nanIdx))
        % Find the stretches
        % get vector that contains the indices for each NaN value within
        % inputStruct.values
        nanIndices = find(nanIdx);
        
        % figure out the first index, from this vector nanIndices, that
        % describes the beginning of a run of consecutive NaN values. The
        % first item in nanIdx will always be the beginning of a run but
        % will not have a diff value of 1 there.
        firstIndicesOfRun = [1, (find(diff(nanIndices) ~= 1)+1)]; 
        % figure out the last index, from this vector nanIndices, that
        % describes the end of a run of consecutive NaN values. The last
        % item of nanIdx will always be the end of a run, btu will also no
        % thave a diff value of 1 there.
        lastIndicesOfRun = [(find(diff(nanIndices) ~= 1)-1), length(nanIndices)]; 
        
        % loop over each run of NaNs
        for ff = 1:length(firstIndicesOfRun)
                
                % identify indices corresponding to the entire runs of
                % NaNs. These indices are relative to inputStruct.
                nanRuns{ff} = nanIndices(firstIndicesOfRun(ff)):nanIndices(lastIndicesOfRun(ff));
                
                % if the length of the NaN run is too long, then censor
                % output after introducing appropriate lag
                if inputStruct.timebase(nanRuns{ff}(end)) - inputStruct.timebase(nanRuns{ff}(1)) > p.Results.durationMsecsOfNansToCensor
                    outputStruct.values(ii,nanRuns{ff}(1)+lagToPeakIndex:nanRuns{ff}(end)+lagToPeakIndex) = NaN;
                end
        end

    end
    
    
end



