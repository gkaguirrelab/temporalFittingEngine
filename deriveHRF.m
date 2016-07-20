function [HRF,fSet,betaValues,DesignMatrix,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs)

% Derives a haemodynamic response function (HRF) using fMRI time-series
% data and an input matrix of event times.
%
%   Usage:
%   [HRF,betaValues,fSet,estTC,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs)
%
%   Inputs:
%   timeSeries      - matrix of time-series data (TR x N)
%   eventTimes      - matrix of event times (msec)
%   sampT           - sampling period (msec) e.g. TR [default = 1000]
%   HRFdur          - duration of HRF to be modeled (msec) [default = 34000]
%   numFreqs        - number of frequencies to be modeled (scalar value)
%
%   note:
%   the value of 'numFreqs' may be reduced when the Fourier set is
%   downsampled to the resolution of the input time-series
%
%   Outputs:
%   HRF             - estimated HRF
%   betaValues      - beta weights for the linearly indepedent covariates
%   fSet            - linearly indepedent Fourier set (msec resolution)
%   DesignMatrix    - linearly indepedent Fourier set (timeSeries resolution)
%   numCov          - number of independent covariates in the Fourier set
%
%   Written by Andrew S Bock Jul 2016

%% Set defaults
if ~exist('HRFdur','var') || isempty(HRFdur)
    HRFdur      = 34000; % msec
end
if ~exist('sampT','var') || isempty(sampT)
    sampT       = 1000; % msec
end
% For Phase Shift
floorTimes = eventTimes ./ 1000 ;               % Floor works with decimals
RndTimes = floor(floorTimes) ;
RndTimes = RndTimes .* 1000 ;
Shift_Amount = eventTimes - RndTimes ;          % Phase Shift Amount
% For Fourier set below
t               = linspace(0,HRFdur-1,HRFdur);  % Create Time (msec) Array
fSet            = zeros(HRFdur,(2*numFreqs)+1); % Create blank Fourier Set (add one for the dc)
%% Create Fourier Set
ct = 1;
fSet(:,ct)      = ones(HRFdur,1);               % Create DC component
for i = 1:numFreqs
    ct = ct + 1;
    fSet(:,ct)  = sin(t/HRFdur*2*pi*i);         % Create Sin waves for each Fq
    ct = ct + 1;
    fSet(:,ct)  = cos(t/HRFdur*2*pi*i);         % Create Cos waves for each Fq
end
% Downsample the Fourier Set
DfSet           = resample(fSet,1,sampT);
% Only keep the linearly independent covariates
[~,goodCovs]    = indMat(DfSet);
fSet            = fSet(:,goodCovs);
fDims           = size(fSet);
numCov          = fDims(2);                     % number of covariates
%% Create the design matrix
msecTC          = size(timeSeries,1)*sampT;     % length of time-series (msec)
tempMatrix      = zeros(msecTC+HRFdur,fDims(2));
for i = 1:length(eventTimes)
    thisBlock   = eventTimes(i)+(0:HRFdur - 1);
    % Phase Shift
    fSet_PS = circshift(fSet,Shift_Amount(i)) ; % Create Phase Shifted Fourier Set
    % DC component -- column of 1's
    tempMatrix(thisBlock,1) = fSet_PS(1:size(fSet_PS,1),1);
    % Add each event (combines [sums +] overlapping Fourier sets)
    tempMatrix(thisBlock,2:end) = tempMatrix(thisBlock,2:end) + ...
        fSet_PS(1:size(fSet_PS,1),2:end);
end
% Crop off Excess Rows (outside time-series)
upMatrix        = tempMatrix(1:msecTC,:);
% Downsample design matrix to resolution of time-series data
DesignMatrix    = resample(upMatrix,1,sampT);
%% Run linear regression
betaValues      = DesignMatrix\timeSeries;
%% Get the estimated hrf
HRF             = fSet_PS * betaValues;
