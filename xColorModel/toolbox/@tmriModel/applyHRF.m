function boldResponse = applyHRF(obj,timebase,neuralResponse,HRF,varargin)
% boldResponse = applyHRF(obj,timebase,neuralResponse,HRF,varargin)
% 
% Apply the HRF (contained in the HRF structure) to the neural response to
% obtain the BOLD response.  For right now, we'll assume that the HRF
% structure contains a single field, hrf, on the same timebase spacing as
% the neural response.
%
% Inputs:
%   timebase - times on which data/model predictions exist (in seconds)
%   neuralResponse - neural response on timebase
%   HRF - structure containing info about the HRF
% 
% Optional key/value pairs

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('timebase',@isnumeric);
p.addRequired('neuralResponse',@isnumeric);
p.addRequired('HRF',@isstruct);
p.parse(timebase,neuralResponse,HRF,varargin{:});

%% Sample at Points t
stimulusUpsampled = interp1(stimTimeSamples,stimulus,t,'linear','extrap') ;
stimulusUpsampled = stimulusUpsampled(1:length(t)) ;

% Convolve Stimulus with HRF to get Regressor
regressorPreCut = conv(stimulusUpsampled,hrf) ;

% Cut off extra Conv values --( Need to look more into this. Conv is
% wierd in Matlab)
regressorUpsampled = regressorPreCut(1:length(stimulusUpsampled)) ;
regressor = interp1(t,regressorUpsampled,stimTimeSamples,'linear','extrap') ;