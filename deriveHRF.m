function deriveHRF(timeSeries,eventTimes,sampT,HRFdur)
% <NEED TO ADD DESCRIPTION ONCE FINALIZED>
%
%   Usage:
%       <NEED TO ADD USAGE ONCE FINALIZED>
%
%   Inputs:
%   timeSeries  - matrix of time-series data (TR x N)
%   eventTimes  - matrix of event times (msec)
%   sampT       - sampling period (msec) e.g. TR
%   HRFdur      - duration of HRF to be modeled (msec)

%% Set defaults
if ~exist('HRFdur','var') || isempty(HRFdur)
    HRFdur      = 33000; % seconds
end
if ~exist('sampT','var') || isempty(sampT)
    sampT       = 1000; % 1 second
end
% For Fourier set below
downHRFdur      = HRFdur / 1000;                % Downsampled HRF (sec)
t               = linspace(0,HRFdur-1,HRFdur);  % Create Time (msec) Array
m               = zeros(downHRFdur,HRFdur);     % Create blank Fourier Set
m(1,:)          = ones(HRFdur,1);               % Create DC component
%% Create Fourier Set
% if the HRF is an odd number of msecs
if mod(HRFdur,2) == 1
    for i = 1:(downHRFdur-1)/2
        m(i*2,:)   = sin(t/HRFdur*2*pi*i);      % Create Sin waves for each Fq
        m(i*2+1,:) = cos(t/HRFdur*2*pi*i);      % Create Cos waves for each Fq
    end
    % if the HRF is an even number of msecs
elseif mod(HRFdur,2) == 0
    for i = 1:floor(downHRFdur/2)-1
        m(i*2,:)   = sin(t/(HRFdur-1)*2*pi*i);  % Create Sin waves for each Fq
        m(i*2+1,:) = cos(t/(HRFdur-1)*2*pi*i);  % Create Cos waves for each Fq
    end
end
% True Nyquist Frequency
m(downHRFdur,:)     = sin(t/(HRFdur-1)*2*pi*(downHRFdur/2));
% Rotate matrix
m               = m';
%% Create the design matrix
msecTC          = size(timeSeries,1)*sampT; % length of time-series (msec) 
tempMatrix      = zeros(msecTC+HRFdur,downHRFdur);
for i = 1:length(eventTimes)
    % DC component -- column of 1's
    tempMatrix(eventTimes(i)+(0:HRFdur - 1),1) = m(1:size(m,1),1);
    % Add each event (combines overlapping Fourier sets)
    tempMatrix(eventTimes(i)+(0:HRFdur - 1),2:downHRFdur) = ...
        tempMatrix(eventTimes(i)+(0:HRFdur - 1),2:downHRFdur) + ...
        m(1:size(m,1),2:downHRFdur);
end
% Crop off Excess Rows (outside time-series)
upMatrix        = tempMatrix(1:msecTC,:);
% Downsample design matrix to resolution of time-series data
DesignMatrix    = downsample(upMatrix,sampT);
%% Run linear regression
b               = DesignMatrix\timeSeries;
%% Get the hrf
hrf = b'*DesignMatrix';
%% plot data
x = 0:size(timeSeries,1)-1;
figure;plot(x,hrf,'r',x,timeSeries,'k',x,timeSeries,'ko');
axis square;
xlabel('Time (TRs)','FontSize',20);
ylabel('Percent Signal Change','FontSize',20);
legend({'Fourier model' 'Simulated data'},'FontSize',20,'Location','NorthEastOutside');