%% MAKES FULL SET OF BDCM PLOTS

% make sure to add mriTemporalFitting and all subdirectories to path. that
% is, cd into mriTemporalFitting, then call addpath(genpath(pwd))

%% LOAD APPROPRIATE FIT RESULTS STRUCT FOR THE SUBJECT OF INTEREST

subjName = 'gka1';

if strcmp(subjName,'gka1')
    load('gka1_fitResults');
elseif strcmp(subjName,'asb1')
    load('asb1_fitResults');
else
   error('makeFullSetOfPlots: pick valid subject'); 
end

% call huge plotting function
fullSetOfPlotsFunc(fitResults);