function regressor = createRegressor(stimulus,stimTimeSamples,hrf,t)

% function regressor = createRegressor(stimulus,stimTimeSamples,hrf,t)
%
% upsamples the stimulus to the same resolution as HRF, then convolves

% Sample at Points t
stimulusUpsampled = interp1(stimTimeSamples,stimulus,t,'linear','extrap') ;
stimulusUpsampled = stimulusUpsampled(1:length(t)) ;

% Convolve Stimulus with HRF to get Regressor
regressorPreCut = conv(stimulusUpsampled,hrf) ;

% Cut off extra Conv values --( Need to look more into this. Conv is
% wierd in Matlab)
regressorUpsampled = regressorPreCut(1:length(stimulusUpsampled)) ;
regressor = interp1(t,regressorUpsampled,stimTimeSamples,'linear','extrap') ;

gribble = 1;