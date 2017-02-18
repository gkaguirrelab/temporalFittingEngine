function [ kernelStructOut ] = prepareHRFKernel( kernelStructIn )
% Modify the passed kernel structure to make it suitable for use as a
% convolution kernel, assuming that the kernel is a BOLD hemodynamic
% response function

% Copy over all fields from input to output
kernelStructOut = kernelStructIn;

% Obtain the deltaT of the kernelStruct timebase
check = diff(kernelStructIn.timebase);
kernelDeltaT = check(1);

% Set the initial value of the kernel to zero
kernelStructOut.values=kernelStructOut.values-kernelStructOut.values(1);

% scale the kernel to preserve amplitude following convolution
kernelStructOut.values=kernelStructOut.values/abs((sum(kernelStructOut.values)*kernelDeltaT));

% Issue a warning that this function is deprecated, and should be replaced
% by setting the initial value of the kernel manually in the calling code,
% and then using normalizeKernelAmplitude() to scale the kernel.
warning('** prepareHRFKernel is deprecated. Replace with call to normalizeKernelAmplitude, and manual setting of initial kernel value');

end

