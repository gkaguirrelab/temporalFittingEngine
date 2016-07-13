function [hrf, timeSeriesNoAttn] = Derive_HRF_using_Fourier(originalTimeSamples,timeSeries,EventStartTimes,HRFduration,HRFsampleInterval)

% infers HRF from Event (Attention Task)
%
% originalTimeSamples: time samples for time series
% timeSeries         : self-explanatory
% attnStartTimes     : when each event began
% HRFduration        : duration of the HRF
% HRFsampleInternal  : TR

% NORMALIZE TIME SAMPLES TO 1
originalTimeSamplesNormalized = originalTimeSamples.*(1./HRFsampleInterval);

%UP-sampled Time Series
UPS_originalTimeSamplesNormalized = ...         % up-sample Time Samples
originalTimeSamplesNormalized * 1000 ;
dt = 1;                                         % (1ms) UPsampled Resolution
UP_R = 1000 ;                                   % (1ms) UPsampled Resolution
UPS_TimeSamples = ...                           % Fill with 1ms Resolution
(0:dt:UPS_originalTimeSamplesNormalized(length(UPS_originalTimeSamplesNormalized))) ;

% Event Time Rounded-- For Phase Shifting
RndEventTimes = floor(EventStartTimes) ;
RndEventTimes = RndEventTimes * UP_R ;

% UP-sample Event TIme
EventStartTimes = round(EventStartTimes,3) ;    % First Round nearest 3rd decimal
EventStartTimes = EventStartTimes * UP_R ;
Shift_Amount = EventStartTimes - RndEventTimes ;

%% Fourier Set

% n = HRF duration-- Time Set (seconds)
% m = Fourier Set Matrix
% t = Time Points (HRF duration)
HRF_TR = HRFduration *UP_R ;                    % Name conventions for waves

% n -- Is odd -------------------------------------------------------------
if mod(HRF_TR,2) == 1                   

    t      = linspace(0,HRF_TR-1,HRF_TR);       % Create Time(s) Array
    m      = zeros(HRFduration,HRF_TR);         % Create blank Fourier Set
    m(1,:) = ones(HRF_TR,1);                    % Create DC component

    for i = 1:(HRFduration-1)/2
        m(i*2,:)   = sin(t/HRF_TR*2*pi*i);      % Create Sin waves for each Fq 
        m(i*2+1,:) = cos(t/HRF_TR*2*pi*i);      % Create Cos waves for each Fq
    end

% n -- Is even ------------------------------------------------------------
elseif mod(HRF_TR,2) == 0                    

    t      = linspace(0,HRF_TR-1,HRF_TR);       % Create Time(s) Array
    m      = zeros(HRFduration,HRF_TR);         % Create blank Fourier Set
    m(1,:) = ones(HRF_TR,1);                    % Create DC component

    for i = 1:floor(HRFduration/2)-1
        m(i*2,:)   = sin(t/(HRF_TR-1)*2*pi*i);  % Create Sin waves for each Fq
        m(i*2+1,:) = cos(t/(HRF_TR-1)*2*pi*i);  % Create Cos waves for each Fq
    end
    
else
% Leave Fourier set Blank -------------------------------------------------
    t = linspace(0,HRF_TR-1,HRF_TR);            % Create Time(s) Array
    m = zeros(HRFduration,HRF_TR);              % Create blank Fourier Set
end

% True Nyquist Frequency --------------------------------------------------
    m(HRFduration,:) = sin(t/(HRF_TR-1)*2*pi*(HRFduration/2));

% ********************************* Leaving IN True Nyquist

% Change Dimensions
m = m' ;
% m = m(:,1:size(m,2)-1);        % In case we want to leave True Nyquis OUT


%% Design Matrix (TR 1ms) -------------------------------------------------

% Create matrix
TimeSeriesMatrix = zeros(length(UPS_TimeSamples)+(HRF_TR),HRFduration) ;

% Create Design Matrix at UP-sampled Event Times & Phase Shift
for i = 1:length(EventStartTimes)
    % DC component -- column of 1's
    TimeSeriesMatrix(EventStartTimes(i)+(0:HRF_TR - 1- Shift_Amount(i)),1) = ...
    m(1:length(m)-Shift_Amount(i),1) ;    

    % Waves at Event Time (ms resolution)
    TimeSeriesMatrix(EventStartTimes(i)+(0:HRF_TR - 1 - Shift_Amount(i)),2:HRFduration-1) = ...
    TimeSeriesMatrix(EventStartTimes(i)+(0:HRF_TR - 1 - Shift_Amount(i)),2:HRFduration-1) + ...
    m(1:length(m)-Shift_Amount(i),2:HRFduration-1) ;

% ---------------- Phase Shift --------------------------------------------

    % Wrap Around to Phase Shift DC Component
    TimeSeriesMatrix(EventStartTimes(i) - Shift_Amount(i) +(0:Shift_Amount(i)-1),1) = ...
    m(length(m)-Shift_Amount(i):length(m)-1,1) ;      

    % Wrap Around to Phase Shift for data resolution
    TimeSeriesMatrix(EventStartTimes(i) - Shift_Amount(i) +(0:Shift_Amount(i)-1),2:HRFduration-1) = ...
    TimeSeriesMatrix(EventStartTimes(i) - Shift_Amount(i) +(0:Shift_Amount(i)-1),2:HRFduration-1) + ...
    m(length(m)-Shift_Amount(i):length(m)-1,2:HRFduration-1) ;   
    
end 

% Crop off Excess Rows-- Leave length of Original Time series
DesignMatrix = TimeSeriesMatrix(1:length(UPS_TimeSamples),:) ;

% Crop off Excess Column-- discarded True Nyquist ***** KEPT NYQUIST
%DesignMatrix = DesignMatrix(:,1:size(DesignMatrix,2)-1) ;
% ^^ This Fixed the Rank Deficiency. NO need to Orthogonalize ^^

%% Downsample to 1s Resolution

% DesignMatrix = downsample(DesignMatrix,UP_R) ;        % Down-sample data Resolution
% Interpolate Design Matrix to downsample
DesignMatrix = interp1(UPS_TimeSamples,DesignMatrix,UPS_originalTimeSamplesNormalized) ;

%% GET BETA VALUES & HRF

% suppress warnings regarding a rank deficient matrix. It's OK 
warning('off','MATLAB:rankDeficientMatrix')
betaValues = DesignMatrix\timeSeries';
warning('on','MATLAB:rankDeficientMatrix')

% HRF -- Fourier Transform (Reconstruct waves with Beta values)
betaValues = betaValues' ;

% Design Matrix with NO Overlaping Event Time
HRF_No_Ovrlap = [] ;

% Down-Sample Fourier 'Cassette' (m - variable)
m_SR = 1:dt:length(m) ;                 % m- Sample Rate
m_SR = m_SR' ;
m_DSR = UP_R:UP_R:length(m) ;           % m- Down Sampling Rate
m_DS = interp1(m_SR,m,m_DSR) ;          % Interpolate Fourier Cassete

% Weight single Fourier Set
for i = 1:length(betaValues)
    HRF_No_Ovrlap(:,i) = m_DS(:,i) .* betaValues(:,i) ;
end

% Reconstruct HRF from single Fourier Set
for j = 1:length(m_DS)
    HRF_No_Ovrlap(j) = sum(HRF_No_Ovrlap(j,:)) ;
end

% Save the first column -- Summed & weighted Fourier Set.
HRF_No_Ovrlap = HRF_No_Ovrlap(:,1) ;    % Take 1st column (Reconstructed hrf)
HRF_No_Ovrlap = HRF_No_Ovrlap' ;
hrf           = HRF_No_Ovrlap ;         % Compatibility with rest of code

%% REGRESS EVENT FROM TIMES SERIES

HRF_Matrix = [] ;
% Weight each wave with its corresponidng beta value
for i = 1:length(betaValues)
    HRF_Matrix(:,i) = DesignMatrix(:,i) .* betaValues(:,i) ;
end

% Reconstruct HRF with weigthed beta values (From every Fourier Set)
HRF = [] ;
for i = 1:length(UPS_originalTimeSamplesNormalized)
    HRF(i) = sum(HRF_Matrix(i,:)) ;
end 

% Timeseries With HRF substracted
timeSeriesNoAttn = timeSeries - HRF ;

Graves = 1;