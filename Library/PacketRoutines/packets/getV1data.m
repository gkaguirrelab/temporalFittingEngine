% Script to load fMRI data from V1 from an example TOME subject

error('This function is deprecated. Define packets within your project, not using the tfe');

%% set defaults
outDir = '~';

%% Load the relevant data
sessionDir      = '/data/jag/TOME/TOME_3001/081916b';
boldDirs        = find_bold(sessionDir); 
f               = load_nifti(fullfile(sessionDir,boldDirs{1},'s5.wdrf.tf.surf.lh.nii.gz'));
a               = load_nifti(fullfile(sessionDir,'anat_templates','lh.areas.anat.nii.gz'));
%% Save the V1 data
V1ind           = find(abs(a.vol)==1);
tmp             = squeeze(f.vol(V1ind,:,:,:));
V1tc            = convert_to_psc(tmp);
save(fullfile(outDir,'V1tc.mat'),'V1tc');