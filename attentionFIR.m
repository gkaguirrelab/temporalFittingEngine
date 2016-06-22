function [hrf, timeSeriesNoAttn] = attentionFIR(originalTimeSamples,timeSeries,attnStartTimes,HRFduration,HRFsampleInterval)

% function [hrf, timeSeriesNoAttn] = attentionFIR(originalTimeSamples,timeSeries,attnStartTimes,HRFduration,HRFsampleInterval)
%
% infers HRF from attention task
%
% originalTimeSamples: time samples for time series
% timeSeries         : self-explanatory
% attnStartTimes     : when each attention task began
% HRFduration        : duration of the HRF
% HRFsampleInternal  : TR

% NORMALIZE TIME SAMPLES TO 1
originalTimeSamplesNormalized = originalTimeSamples.*(1./HRFsampleInterval);

% HELPER VECTOR FOR SHIFT
timeShiftScale = [1:HRFduration./HRFsampleInterval];

% MATRIX OF START TIMES TO BE SHIFTED
timesToBeShifted = repmat(attnStartTimes,[1 HRFduration]);

% IMPLEMENT SHIFT
timeShiftedMatrix = timesToBeShifted+repmat(timeShiftScale,[size(timesToBeShifted,1) 1]).*HRFsampleInterval;

timeShiftedMatrix = round(timeShiftedMatrix);

% CREATES DESIGN MATRIX
for i = 1:length(timeShiftScale);
   impulseVector = double(ismember(originalTimeSamplesNormalized,timeShiftedMatrix(:,i)));
%   designMatrix(:,i) = impulseVector;
   designMatrix(:,i) = impulseVector - mean(impulseVector);
end

% ADD THE REGRESSOR FOR THE MEAN
designMatrix = [ones([size(designMatrix,1) 1]) designMatrix];

% GET BETA VALUES
betaValues = designMatrix\timeSeries';

% HRF IS ALL BETA VALUES EXCEPT FOR THE MEAN VECTOR
hrf = betaValues(2:length(betaValues))';

% REGRESS ATTENTION FROM TIMES SERIES
timeSeriesNoAttn = timeSeries - sum(designMatrix(:,2:size(designMatrix,2)).*repmat(hrf,[size(designMatrix,1) 1]),2)';

gribble = 1;