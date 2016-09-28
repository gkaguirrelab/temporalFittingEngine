function tfeConvolutionDemo
% tfeConvolutionDemo
%
% Demonstrate and test the applyKernel method of the TFE object
%
% 6/26/16  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Construct a model object
%
% We can't instantiate the parent class, so we pick one of the
% subclasses. It shouldn't matter which one, as we are not going
% to use an
tfe = tfeQCM('verbosity','none');

%% Create a response structure to convolve with
%
% We start with something very simple, so it is easy
% to check that the convolution is working right.
responseDeltaT = 1;
totalResponseTime = 100;
responseStruct.timebase = 0:responseDeltaT:totalResponseTime;
responseStruct.values = zeros(size(responseStruct.timebase));
responseNonZeroEntries = [1 round(totalResponseTime/2) totalResponseTime+1];
responseStruct.values(responseNonZeroEntries) = 1;

%% Create a kernel struct
%
% Again start with something simple
kernelDeltaT = 1;
totalKernelTime = 100;
kernelStruct.timebase = 0:kernelDeltaT:totalKernelTime;
kernelStruct.values = zeros(size(kernelStruct.timebase));
kernelNonZeroEntries = [1];
kernelStruct.values(kernelNonZeroEntries) = 1;

%% Apply kernel to response
convResponseStruct = tfe.applyKernel(responseStruct,kernelStruct);

%% Make a plot
figure; clf;
subplot(3,1,1); hold on
plot(responseStruct.timebase,responseStruct.values,'ro','MarkerFaceColor','r','MarkerSize',8);
subplot(3,1,2); hold on
plot(kernelStruct.timebase,kernelStruct.values,'go','MarkerFaceColor','g','MarkerSize',8);
subplot(3,1,3); hold on
plot(convResponseStruct.timebase,convResponseStruct.values,'bo','MarkerFaceColor','b','MarkerSize',8);

end

