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

% n -- Is odd -------------------------------------------------------------
if mod(n,2) == 1                   

    t      = linspace(0,n-1,n);             % Create Time(s) Array
    m      = zeros(n,n);                    % Create blank Fourier Set
    m(1,:) = ones(n,1);                     % Create DC component

    for i = 1:floor(n/2)-1
        m(i*2,:)   = sin(t/n*2*pi*i);       % Create Sin waves for each Fq 
        m(i*2+1,:) = cos(t/n*2*pi*i);       % Create Cos waves for each Fq
    end

% n -- Is even ------------------------------------------------------------
elseif mod(n,2) == 0                    

    t      = linspace(0,n-1,n);             % Create Time(s) Array
    m      = zeros(n,n);                    % Create blank Fourier Set
    m(1,:) = ones(n,1);                     % Create DC component

    for i = 1:floor(n/2)-1
        m(i*2,:)   = sin(t/(n-1)*2*pi*i);   % Create Sin waves for each Fq
        m(i*2+1,:) = cos(t/(n-1)*2*pi*i);   % Create Cos waves for each Fq
    end
    
else
% Leave Fourier set Blank -------------------------------------------------
    t = linspace(0,n-1,n);                  % Create Time(s) Array
    m = zeros(n,n);                         % Create blank Fourier Set
end

% Nyquist Frequency -------------------------------------------------------
m(n,:) = sin(t/(n-1)*2*pi*(n/2));

% Change Dimensions
m = m' ;

%% Design Matrix

TimeSeriesMatrix = zeros(length(originalTimeSamplesNormalized),HRFduration) ;

for i = 1:length(RndAttnTimes)
    % DC component -- column of 1's
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),1) = m(1:length(m),1) ;    

    % Everything else
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),2:HRFduration) =  ...
    TimeSeriesMatrix(RndAttnTimes(i)+(0:HRFduration - 1),2:HRFduration) + m(1:length(m),2:HRFduration) ;
    
    % Crop off Extra points. Leave length of Original Time series
    DesignMatrix = TimeSeriesMatrix(1:length(originalTimeSamplesNormalized),:) ;
end 


%% GET BETA VALUES & HRF
betaValues = DesignMatrix\timeSeries';

% HRF -- Fourier Transform (Reconstruct waves with Beta values)
HRF_Matrix = [] ;
betaValues = betaValues' ;

% Weight each wave with its corresponidng beta value
for i = 1:length(betaValues)
    HRF_Matrix(:,i) = DesignMatrix(:,i) .* betaValues(:,i) ;
end

% Reconstruct HRF with weigthed beta values (From every Fourier Set)
HRF = [] ;
for i = 1:length(originalTimeSamplesNormalized)
    HRF(i) = sum(HRF_Matrix(i,:)) ;
end 

%% Design Matrix with NO Overlaping attention tasks
HRF_No_Ovrlap = [] ;

% Weight single Fourier Set
for i = 1:length(betaValues)
    HRF_No_Ovrlap(:,i) = m(:,i) .* betaValues(:,i) ;
end

% Reconstruct HRF from single Fourier Set
for j = 1:length(m)
    HRF_No_Ovrlap(j) = sum(HRF_No_Ovrlap(j,:)) ;
end
% Save the first column -- Summed & weighted Fourier Set.
HRF_No_Ovrlap = HRF_No_Ovrlap(:,1) ;
HRF_No_Ovrlap = HRF_No_Ovrlap' ;

%% REGRESS ATTENTION FROM TIMES SERIES
timeSeriesNoAttn = timeSeries - sum(DesignMatrix(:,2:size(DesignMatrix,2)).*repmat(HRF_No_Ovrlap(:,2:size(HRF_No_Ovrlap,2)),[size(DesignMatrix,1) 1]),2)';

gribble = 1;