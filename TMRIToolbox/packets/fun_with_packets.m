% Example script outlining how to use 'makePacket'
%
%   Written by Andrew S Bock Sep 2016

%% Set initial params
params.packetType       = 'bold';
params.sessionDir       = '/data/jag/MELA/MelanopsinMR/HERO_asb1/032416';
params.runNum           = 1;
params.stimulusFile     = fullfile(params.sessionDir,'MatFiles/HERO_asb1-MelanopsinMRMaxMel-01.mat');
params.responseFile     = fullfile(params.sessionDir,'Series_012_fMRI_MaxMelPulse_A_AP_run01/wdrf.tf.nii.gz');
%% load the response file
tmp                     = load_nifti(params.responseFile);
TR                      = tmp.pixdim(5)/1000;
runDur                  = size(tmp.vol,4);
params.respTimeBase     = 0:TR:(runDur*TR)-TR;
%% If 'bold', get HRF
if strcmp(params.packetType,'bold')
    params.hrfFile      = fullfile(params.sessionDir,'HRF','V1.mat');
end
%% Run a 'loop'
params.timeSeries       = squeeze(tmp.vol(50,50,30,:))';
packet                  = makePacket(params);
eventNum                = 1; % first stimulus event
[B,R2]                  = dummyFit(packet,eventNum);


