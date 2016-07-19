function [HRF,betaValues,fSet,estTC,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs)
% Derives a haemodynamic response function (HRF) using fMRI time-series
% data and an input matrix of event times. 
%
%   Usage:
%       [HRF,betaValues,fSet,estTC,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs)
%
%   Inputs:
%   timeSeries  - matrix of time-series data (TR x N)
%   eventTimes  - matrix of event times (msec)
%   sampT       - sampling period (msec) e.g. TR [default = 1000]
%   HRFdur      - duration of HRF to be modeled (msec) [default = 34000]
%   numFreqs    - number of frequencies to be modeled
%

%% Set defaults
if ~exist('HRFdur','var') || isempty(HRFdur)
    HRFdur      = 34000; % msec
end
if ~exist('sampT','var') || isempty(sampT)
    sampT       = 1000; % msec
end
% For Fourier set below
t               = linspace(0,HRFdur-1,HRFdur);  % Create Time (msec) Array
fSet            = zeros(HRFdur,(2*numFreqs)+1); % Create blank Fourier Set (add one for the dc)
%% Create Fourier Set
ct = 1;
fSet(:,ct)       = ones(HRFdur,1);               % Create DC component
for i = 1:numFreqs
    ct = ct + 1;
    fSet(:,ct) = sin(t/HRFdur*2*pi*i);      % Create Sin waves for each Fq
    ct = ct + 1;
    fSet(:,ct) = cos(t/HRFdur*2*pi*i);      % Create Cos waves for each Fq
end
fDims = size(fSet);
numCov = (fDims(2) - 1) / 2;
%% Plot Fourier set
% figure;
% for i = 1:size(fSet,2);
%     subplot(7,10,i);
%     plot(fSet(:,i));
%     xlim([0 33000]);
% end
%% Create the design matrix
msecTC          = size(timeSeries,1)*sampT; % length of time-series (msec)
tempMatrix      = zeros(msecTC+HRFdur,fDims(2)); 
for i = 1:length(eventTimes)
    % DC component -- column of 1's
    tempMatrix(eventTimes(i)+(0:HRFdur - 1),1) = fSet(1:size(fSet,1),1);
    % Add each event (combines overlapping Fourier sets)
    tempMatrix(eventTimes(i)+(0:HRFdur - 1),2:end) = ...
        tempMatrix(eventTimes(i)+(0:HRFdur - 1),2:end) + ...
        fSet(1:size(fSet,1),2:end);
end
% Crop off Excess Rows (outside time-series)
upMatrix        = tempMatrix(1:msecTC,:);
% Downsample design matrix to resolution of time-series data
DesignMatrix    = downsample(upMatrix,sampT);
%% Run linear regression
betaValues      = DesignMatrix\timeSeries;
%% Get the estimated time-series
estTC = betaValues'*DesignMatrix';
%%
%% Get the estimated hrf
HRF = fSet * betaValues; % ignore dc beta
return;
figure;plot(HRF);
%% plot data
x = 0:size(timeSeries,1)-1;
figure;plot(x,estTC,'r',x,timeSeries,'k',x,timeSeries,'ko');
axis square;
xlabel('Time (TRs)','FontSize',20);
ylabel('Percent Signal Change','FontSize',20);
legend({'Fourier model' 'Simulated data'},'FontSize',20,'Location','NorthEastOutside');