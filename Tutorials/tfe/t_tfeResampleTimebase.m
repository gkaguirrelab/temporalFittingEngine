function validationData = t_tfeResampleTimebase(varargin)
% validationData = t_tfeResampleTimebase(varargin)
%
% Assessments of the resampleTimebase method in the parent tfe class
%
% Optional key/value pairs
%  'originalDeltaT' - value (default 1).  Response temporal sampling in msec.
%  'originalDuration' - value (default 100).  Response duration in msec.
%  'newDeltaT' - value (default 4).  New timebase temporal sampling in msec.
%  'generatePlots' - true/fale (default true).  Make plots?
%  'method' - Resampling method that is passed on to resampleTimebase

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addParameter('originalDeltaT',1,@isnumeric);
p.addParameter('originalDuration',100,@isnumeric);
p.addParameter('newDeltaT',4,@isnumeric);
p.addParameter('newDuration',100,@isnumeric);
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});

%% Construct the model object
%
% We have to pick some subclass to define the model, although this
% shouldn't matter for the validation tests here
tfe = tfeIAMP('verbosity','high');

%% Create a original structure to change the timebase of
originalStruct.timebase = 0;
while (originalStruct.timebase(end) < p.Results.originalDuration)
    originalStruct.timebase(end+1) = originalStruct.timebase(end)+p.Results.originalDeltaT;
end
originalStruct.values = zeros(size(originalStruct.timebase));
originalNonZeroEntries = [1 round(p.Results.originalDuration/4) round(p.Results.originalDuration/4)+1 round(p.Results.originalDuration/4)+2 round(p.Results.originalDuration/2) p.Results.originalDuration+1];
originalStruct.values(originalNonZeroEntries) = 1;

%% Create a new timebase
newTimebase = 0;
while (newTimebase(end) < p.Results.newDuration)
    newTimebase(end+1) = newTimebase(end)+p.Results.newDeltaT;
end

%% Resample
[resampledStruct] = tfe.resampleTimebase(originalStruct,newTimebase,varargin{:});

%% Make a plot
if (p.Results.generatePlots)
    figure; clf;
    subplot(2,1,1); hold on
    plot(originalStruct.timebase,originalStruct.values,'ro','MarkerFaceColor','r','MarkerSize',8);
    xlabel('Time (msec)'); ylabel('Response'); title('Original');

    subplot(2,1,2); hold on
    plot(resampledStruct.timebase,resampledStruct.values,'bo','MarkerFaceColor','b','MarkerSize',8);
    xlabel('Time (msec)'); ylabel('Response');  title('Resampled');
end

%% Set validation data for return
if (nargout > 0)
    validationData.originalStruct = originalStruct;
    validationData.resampledKernelStruct = resampledStruct;
end

end

