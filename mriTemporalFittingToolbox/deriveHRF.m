function [HRF,cleanData,fSet,betaValues,DesignMatrix,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs,modelType,phaseShift)

% Derives a haemodynamic response function (HRF) using fMRI time-series
% data and an input matrix of event times.
%
%   Usage:
%   [HRF,cleanData,fSet,betaValues,DesignMatrix,numCov] = ...
%       deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs,modelType,phaseShift)
%
%   Inputs:
%   timeSeries      - matrix of time-series data (TR x N)
%   eventTimes      - matrix of event times (msec)
%   sampT           - sampling period (msec) e.g. TR [default = 1000]
%   HRFdur          - duration of HRF to be modeled (msec) [default = 32000]
%   numFreqs        - number of frequencies to be modeled (scalar value)
%
%   note:
%   the value of 'numFreqs' may be reduced when the Fourier set is
%   downsampled to the resolution of the input time-series
%
%   Outputs:
%   HRF             - estimated HRF
%   cleanData       - input timeseries with attention event responses regressed out
%   betaValues      - beta weights for the linearly indepedent covariates
%   fSet            - linearly indepedent Fourier set (msec resolution)
%   DesignMatrix    - linearly indepedent Fourier set (timeSeries resolution)
%   numCov          - number of independent covariates in the Fourier set
%
%   Written by Andrew S Bock & Jovan Ortiz Bernhardt - Jul 2016
%   Order of outputs updated by Andrew S Bock Aug 2016

%% Set defaults
if ~exist('sampT','var') || isempty(sampT)
    sampT           = 1000; % msec
end
if ~exist('HRFdur','var') || isempty(HRFdur)
    HRFdur          = 16000; % msec
end
if ~exist('modelType','var') || isempty(modelType)
    modelType       = 'Fourier';
end
if ~exist('phaseShift','var') || isempty(phaseShift)
    phaseShift      = 0;
end
% For Fourier set below
t                   = linspace(0,HRFdur-1,HRFdur);  % Create Time (msec) Array
fSet                = zeros(HRFdur,(2*numFreqs)+1); % Create blank Fourier Set (add one for the dc)
%% Create Fourier Set
ct = 1;
fSet(:,ct)      = ones(HRFdur,1);               % Create DC component
for i = 1:numFreqs
    ct = ct + 1;
    fSet(:,ct)  = sin(t/HRFdur*2*pi*i);         % Create Sin waves for each Fq
    ct = ct + 1;
    fSet(:,ct)  = cos(t/HRFdur*2*pi*i);         % Create Cos waves for each Fq
end
% % Downsample the Fourier Set
DfSet           = resample(fSet,1,sampT);
% Only keep the linearly independent covariates
[~,goodCovs]    = indMat(DfSet);
fSet            = fSet(:,goodCovs);

fDims           = size(fSet);
%% CREATE FIR SET
% initialize as zeros
FIRset = zeros(HRFdur,HRFdur/sampT);
% place spikes at the beginning of event, then at sampT intervals
deltas = (0:sampT:HRFdur)+1; % add 1 to eventTimes (sec -> index)
% loop over indices, place a 1 at each index
for j = 1:length(deltas)-1
    FIRset(deltas(j):deltas(j+1)-1,j) = 1;
end
%% Create the design matrix
msecTC          = size(timeSeries,1)*sampT;     % length of time-series (msec)
tempMatrix      = zeros(msecTC+HRFdur,fDims(2));
for i = 1:length(eventTimes)
    thisBlock   = (eventTimes(i)+1) + (0:HRFdur - 1); % add 1 to eventTimes (sec -> index)
    switch modelType
        case 'Fourier'
            if phaseShift
                % get the shift amount
                FlrTimes = eventTimes ./ 1000 ;
                RndTimes = floor(FlrTimes) ;
                RndTimes = RndTimes .* 1000 ;
                Shift_Amount = eventTimes - RndTimes ;
                % shift the Fourier set
                tmpSet = circshift(fSet,Shift_Amount(i)) ;
                % Add each event (combines overlapping Fourier sets)
                tempMatrix(thisBlock,:) = tempMatrix(thisBlock,:) + ...
                    tmpSet(1:size(tmpSet,1),:);
            else
                % Add each event (combines overlapping Fourier sets)
                tempMatrix(thisBlock,:) = tempMatrix(thisBlock,:) + ...
                    fSet(1:size(fSet,1),:);
            end
        case 'FIR'
            % place the FIR set
            tempMatrix(thisBlock,:) = tempMatrix(thisBlock,:) + FIRset;
    end
end
% Crop off Excess Rows (outside time-series)
upMatrix                = tempMatrix(1:msecTC,:);
% Downsample design matrix to resolution of time-series data
DesignMatrix            = resample(upMatrix,1,sampT);
switch modelType
    case 'Fourier'
        % Only keep the linearly independent covariates
%         [~,goodCovs]    = indMat(DesignMatrix);
%         DesignMatrix    = DesignMatrix(:,goodCovs);
%         fSet            = fSet(:,goodCovs);
    case 'FIR'
        % Mean center the design matrix
        DesignMatrix    = DesignMatrix - repmat(mean(DesignMatrix),size(DesignMatrix,1),1);
end
numCov          = size(DesignMatrix,2); % number of covariates
%% Run linear regression
if strcmp(modelType,'FIR')
    DesignMatrix        = [ones([size(DesignMatrix,1) 1]) DesignMatrix];
end
betaValues              = DesignMatrix\timeSeries;
%% Get the estimated hrf
switch modelType
    case 'Fourier'
        tmp             = fSet * betaValues;
        HRF             = tmp - mean(tmp);
        cleanData       = timeSeries - DesignMatrix*betaValues;
    case 'FIR'
        tmp             = DesignMatrix(:,2:end);
        maxDelta        = max(tmp(:));
        HRF             = betaValues(2:end)*maxDelta - mean(betaValues(2:end)*maxDelta);
        cleanData       = timeSeries - DesignMatrix*betaValues;
end