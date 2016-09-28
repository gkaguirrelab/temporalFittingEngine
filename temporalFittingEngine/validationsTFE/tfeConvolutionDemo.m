function tfeConvolutionDemo(varargin)
% tfeConvolutionDemo
%
% Demonstrate and test the applyKernel method of the TFE object
%
% Optional key/value pairs
%  'responseDeltaT' - value (default 1).  Response temporal sampling in msec.
%  'responseDuration' - value (default 100).  Response duration in msec.
%  'kernelDeltaT' - value (default 1).  Kernel temporal sampling in msec.
%  'kernelDuration' - value (default 100).  Kernel duration in msec.

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('responseDeltaT',1,@isnumeric);
p.addParameter('responseDuration',100,@isnumeric);
p.addParameter('kernelDeltaT',1,@isnumeric);
p.addParameter('kernelDuration',10,@isnumeric);
p.parse(varargin{:});
kernelDeltaT = p.Results.kernelDeltaT;
kernelDuration = p.Results.kernelDuration;
responseDeltaT = p.Results.responseDeltaT;
responseDuration = p.Results.responseDuration;
kernelNonZeroTime = 2;
kernelValue = 2;

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
responseStruct.timebase = 0:responseDeltaT:responseDuration;
if (responseStruct.timebase(end) < responseDuration)
    responseStruct.timebase(end+1) = responseStruct.timebase(end)+responseDeltaT;
end
responseStruct.values = zeros(size(responseStruct.timebase));
responseNonZeroEntries = [1 round(responseDuration/4) round(responseDuration/4)+1 round(responseDuration/4)+2 round(responseDuration/2) responseDuration+1];
responseStruct.values(responseNonZeroEntries) = 1;

%% Create a kernel struct
%
% Again start with something simple
kernelStruct.timebase = 0:kernelDeltaT:kernelDuration;
if (kernelStruct.timebase(end) < kernelDuration)
    kernelStruct.timebase(end+1) = kernelStruct.timebase(end)+kernelDeltaT;
end
kernelStruct.values = zeros(size(kernelStruct.timebase));
index = find(kernelStruct.timebase <= kernelNonZeroTime);
kernelStruct.values(index) = kernelValue;

%% Apply kernel to response
[convResponseStruct,resampledKernelStruct] = tfe.applyKernel(responseStruct,kernelStruct);

%% Make a plot
figure; clf;
subplot(4,1,1); hold on
plot(responseStruct.timebase,responseStruct.values,'ro','MarkerFaceColor','r','MarkerSize',8);
xlabel('Time (msec)'); ylabel('Response'); title('Input');
subplot(4,1,2); hold on
plot(kernelStruct.timebase,kernelStruct.values,'go','MarkerFaceColor','g','MarkerSize',8);
xlabel('Time (msec)'); ylabel('Response');  title('Kernel');
subplot(4,1,3); hold on
plot(resampledKernelStruct.timebase,resampledKernelStruct.values,'ko','MarkerFaceColor','k','MarkerSize',8);
xlabel('Time (msec)'); ylabel('Response');  title('Resampled Kernel');
subplot(4,1,4); hold on
plot(convResponseStruct.timebase,convResponseStruct.values,'bo','MarkerFaceColor','b','MarkerSize',8);
xlabel('Time (msec)'); ylabel('Response');  title('Output');

end

