function [ kernelStructOut ] = normalizeKernelAmplitude( kernelStructIn )
% Modify the passed kernel structure so that the kernel will preserve
% signal amplitude following convolution

% Copy over all fields from input to output
kernelStructOut = kernelStructIn;

% Obtain the deltaT of the kernelStruct timebase
check = diff(kernelStructIn.timebase);
kernelDeltaT = check(1);

% scale the kernel to preserve area following convolution
kernelStructOut.values=kernelStructOut.values/abs((sum(kernelStructOut.values)*kernelDeltaT));

warning('This function is normalizing kernel area, which is usually what is desired. Switch to using the normalizeKernelArea function');

end

