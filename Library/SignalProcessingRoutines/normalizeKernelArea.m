function [ kernelStructOut ] = normalizeKernelArea( kernelStructIn )
% Modify the passed kernel structure so that the kernel will preserve
% signal area following convolution

% Copy over all fields from input to output
kernelStructOut = kernelStructIn;

% Obtain the deltaT of the kernelStruct timebase
check = diff(kernelStructIn.timebase);
kernelDeltaT = check(1);

% scale the kernel to preserve area following convolution
kernelStructOut.values=kernelStructOut.values/(sum(abs(kernelStructOut.values)*kernelDeltaT));

end

