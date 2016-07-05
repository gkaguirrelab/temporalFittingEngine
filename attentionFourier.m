function [hrf, timeSeriesNoAttn] = attentionFourier(originalTimeSamples,timeSeries,attnStartTimes,HRFduration,HRFsampleInterval)

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

% Time Matrix
RndAttnTimes = floor(attnStartTimes) ;

% HELPER VECTOR FOR SHIFT
timeShiftScale = [1:HRFduration./HRFsampleInterval];

%% Fourier Set

% n = HRF duration-- Time Set (seconds)
% m = Fourier Set Matrix
% t = Time Points (HRF duration)
n = HRFduration ;

if mod(n,2) == 1                        % n is odd

    t      = linspace(0,n-1,n);
    m      = zeros(n,n);
    m(1,:) = ones(n,1);

    for i = 1:(n-1)/2
        m(i*2,:)   = sin(t/n*2*pi*i); 
        m(i*2+1,:) = cos(t/n*2*pi*i);
    end
    
elseif mod(n,2) == 0                    % n is even

    t      = linspace(0,n-1,n);
    m      = zeros(n,n);
    m(1,:) = ones(n,1);

    for i = 1:n/2-1
        m(i*2,:)   = sin(t/(n-1)*2*pi*i);
        m(i*2+1,:) = cos(t/(n-1)*2*pi*i);
    end
    
else
    t = linspace(0,n-1,n);
    m = zeros(n,n);
end

m(n,:) = sin(t/(n-1)*2*pi*(n/2));
m = m' ;

%% Design Matrix

TimeSeriesMatrix = zeros(length(originalTimeSamplesNormalized),HRFduration) ;

for i = 1:length(RndAttnTimes)
    % DC component -- column of 1's
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),1) = m(1:length(m),1) ;    

    % Everything else
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),2:HRFduration) =  ...
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),2:HRFduration) + m(1:length(m),2:HRFduration) ;
    
    % Crop off Extra points. Leave at length of Original Time series
    DesignMatrix = TimeSeriesMatrix(1:length(originalTimeSamplesNormalized),:) ;
end 


%% GET BETA VALUES
betaValues = DesignMatrix\timeSeries';

% HRF IS ALL BETA VALUES EXCEPT FOR THE MEAN VECTOR
hrf = betaValues(2:length(betaValues))';

% REGRESS ATTENTION FROM TIMES SERIES
timeSeriesNoAttn = timeSeries - sum(DesignMatrix(:,2:size(DesignMatrix,2)).*repmat(hrf,[size(DesignMatrix,1) 1]),2)';

gribble = 1;