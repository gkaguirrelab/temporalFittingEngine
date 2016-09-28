function outputStruct = applyKernel(obj,inputStruct,kernelStruct,varargin)
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

% Check how many rows are in the the inputStruct values
nRows=size(inputStruct.values,1);

% Check how many timepoints are in the inputStruct timebase
nTimepoints=size(inputStruct.timebase,2);

%% Convolve
%
% But if empty matrix is passed for kernel, return
if (isempty(kernelStruct) || isempty(kernelStruct.values))
    return
end



% Align kernel with 0, if not already done
    zeroAlignedKernel = kernelStruct.values-kernelStruct.values(1);
    
    % The kernel and modelResponseStruct need to be same length: pad with 0's
    convKernel = zeros(1,nTimepoints);
    convKernel(1:length(zeroAlignedKernel)) = zeroAlignedKernel;
 
    % loop over rows for the convolution
for ii=1:nRows
    % Convolve a row of inputStruct.values with the kernel
    valuesRowConv = conv(inputStruct.values(ii,:),convKernel) ;
        % Cut off extra conv values. Need to double check that this is right
    outputStruct.values(ii,:) = valuesRowConv(1:nTimepoints);  
end % loop over rows

end

    
    