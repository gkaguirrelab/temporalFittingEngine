% tests the predictions for the retinotopy model

%% set defaults
dataDir = '/Users/Shared/Matlab/gkaguirrelab/mriTemporalFitting/TMRIToolbox/pRF';

%% load the stimulus file
tmp = load(fullfile(dataDir,'pRFimages.mat'));

%% Binarize the stimulus
stim = tmp.imagesFull;
stim(stim ~=128) = 1;
stim(stim == 128 ) = 0;

%% Make simulated timecourse 
TC = makePrarmsTC(stim,540,540,10);

%% test the fit
R = @(paramsVec)retCorr(stim,TC,paramsVec);
fmincon(R,[500, 500, 8],[],[],[],[],[0 0 0],[1080*1.5,1080*1.5,30])
