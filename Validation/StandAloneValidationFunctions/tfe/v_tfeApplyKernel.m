function validationData = tfeConvolutionDemo(varargin)
% validationData = tfeConvolutionDemo(varargin)
%
% Demonstrate and test the applyKernel method of the TFE object.  This is
% set up so that it may be run as a standalone or incorporated into a
% UnitTestToolbox validation script.
%
% Optional key/value pairs
%  'responseDeltaT' - value (default 1).  Response temporal sampling in msec.
%  'responseDuration' - value (default 100).  Response duration in msec.
%  'kernelDeltaT' - value (default 1).  Kernel temporal sampling in msec.
%  'kernelDuration' - value (default 100).  Kernel duration in msec.
%  'kernelNonZeroTime' - value (default 2). Time over which kernel is
%    non-zero.
%  'kernelValue' - value (default 2).  Non-zero value of kernel.
%  'makePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('responseDeltaT',1,@isnumeric);
p.addParameter('responseDuration',100,@isnumeric);
p.addParameter('kernelDeltaT',1,@isnumeric);
p.addParameter('kernelDuration',10,@isnumeric);
p.addParameter('kernelNonZeroTime',2,@isnumeric);
p.addParameter('kernelValue',2,@isnumeric);
p.addParameter('makePlots',true,@islogical);
p.parse(varargin{:});

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
%
% For reasons I don't understand, with non-integer deltaT there
% is an accumulated spacing error, which is why timebase is pushed
% back to 8 significant places.  This prevents an error when applyKernel
% checks that its inputs are on a regularly spaced timebase.
placesKeep = 4;
responseStruct.timebase = 0;
while (responseStruct.timebase(end) < p.Results.responseDuration)
    responseStruct.timebase(end+1) = responseStruct.timebase(end)+p.Results.responseDeltaT;
end
responseStruct.values = zeros(size(responseStruct.timebase));
responseNonZeroEntries = [1 round(p.Results.responseDuration/4) round(p.Results.responseDuration/4)+1 round(p.Results.responseDuration/4)+2 round(p.Results.responseDuration/2) p.Results.responseDuration+1];
responseStruct.values(responseNonZeroEntries) = 1;

%% Create a kernel struct
%
% Again start with something simple
kernelStruct.timebase = 0;
while (kernelStruct.timebase(end) < p.Results.kernelDuration)
    kernelStruct.timebase(end+1) = kernelStruct.timebase(end)+p.Results.kernelDeltaT;
end
kernelStruct.values = zeros(size(kernelStruct.timebase));
index = find(kernelStruct.timebase <= p.Results.kernelNonZeroTime);
kernelStruct.values(index) = p.Results.kernelValue;

%% Apply kernel to response
[convResponseStruct,resampledKernelStruct] = tfe.applyKernel(responseStruct,kernelStruct);

%% Make a plot
if (p.Results.makePlots)
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
    
    validationData.responseStruct = responseStruct;
    validationData.kernelStruct = kernelStruct;
    validationData.convResponseStruct = convResponseStruct;
    validationData.resampledKernelStruct = resampledKernelStruct;
end

end

