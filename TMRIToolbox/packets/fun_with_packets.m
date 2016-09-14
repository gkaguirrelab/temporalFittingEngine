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
resp                    = load_nifti(params.responseFile);
TR                      = resp.pixdim(5)/1000;
runDur                  = size(resp.vol,4);
params.respTimeBase     = 0:TR:(runDur*TR)-TR;
%% If 'bold', get HRF
if strcmp(params.packetType,'bold')
    params.hrfFile      = fullfile(params.sessionDir,'HRF','V1.mat');
end
%% Loop through the reponse data
volDims                 = size(resp.vol);
flatVol                 = reshape(resp.vol,volDims(1)*volDims(2)*volDims(3),volDims(4));
% Convert to percent signal change
pscVol                  = convert_to_psc(flatVol);
B                       = nan(1,size(pscVol,1));
R2                      = nan(1,size(pscVol,1));
progBar                 = ProgressBar(size(pscVol,1),'looping..');
for i = 1:size(pscVol,1)
    params.timeSeries       = pscVol(i,:);
    packet                  = makePacket(params);
    eventNum                = 1; % first stimulus event
    [B(i),R2(i)]            = dummyFit(packet,eventNum);
    progBar(i);
end


