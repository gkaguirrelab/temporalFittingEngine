%% Notes (25 July)

% add to path
addpath(genpath('/Applications/freesurfer/matlab/'));
addpath(genpath('/Users/Shared/Matlab/gkaguirrelab'));
%addpath(genpath('/Users/jovanortiz/MELA'));

%% Define various variables
% Switches
volType       = ('volume') ;       % surface   | volume

% Variable & paths
session_dir    = '/Users/jovanortiz/MELA/MelanopsinMR/HERO_asb1/040716/' ;
subject_name   = 'HERO_asb1_MaxMel' ;
stimDirs       = listdir(fullfile(session_dir,'Stimuli'),'dirs') ;
stimDir        = fullfile(session_dir,'Stimuli') ;
eventTimesFile = '*MirrorsOffMaxLMS_0Pct_delta_00Sec*' ; % -- Attn Times
%eventTimesFle = '*LMSDirectedSuperMaxLMS_400Pct_delta_00Sec*' ;

HRFdur        = 16000;
numFreqs      = 16;

b              = find_bold(session_dir) ;
runNum         = 1 ;
hemis          = {'lh' 'rh' 'mh'} ;
func           = 'wdrf.tf' ;
tempName       = 'areas.anat.nii.gz'; 

%% Project templates to functional volumes
% Only used once -- (Dont project templates each time)
% Comment & Uncomment next steps accordingly

% interp = 'nearest';
% inv    = 1; % project from anatomical to functional space
% for hh = 1:length(hemis)
%     for i = 1:length(b)
%         invol   = fullfile(session_dir,'anat_templates',[hemis{hh} '.areas.anat.vol.nii.gz']);
%         outvol  = fullfile(session_dir, b{i},[hemis{hh} '.areas.func.vol.nii.gz']) ;
%         targvol = fullfile(session_dir, b{i},[func '.nii.gz']);
%         reg     = fullfile(session_dir, b{i},'func_bbreg.dat');
%         mri_vol2vol(targvol,invol,outvol,reg,interp,inv); % switch invol and targvol, since we're using 'inv';
%     end
% end

%% Load volumes (volume or surface)
%hrf = nan(length(b),16) ;
for ii = 1:length(b)
hemiTCs = [] ;
hh = 3; % choose hemi (1=lh, 2=rh, 3=mh);

% Switch between volume and surface
switch volType
    case 'volume'
        thisMask = fullfile(session_dir, b{ii},[hemis{hh} '.areas.func.vol.nii.gz']);
        thisVol  = fullfile(session_dir, b{ii}, [func '.nii.gz']);
        % Obtain slice timing information
        thisTiming   = load(fullfile(session_dir, b{ii}, 'slicetiming'));
        uniqueTiming = unique(thisTiming);
        
    case 'surface'
        thisMask = fullfile(session_dir, [hemis{hh} '.areas.anat.nii.gz']);
        thisVol  = fullfile(session_dir, b{ii}, [func '.surf.' hemis{hh} '.nii.gz']);
end % End switch

disp(['loading ' thisMask]);            % Display Loading Mask
Mask     = load_nifti(thisMask) ;       % Load Masl
disp(['loading ' thisVol]);             % Display Loading Volume
Volume   = load_nifti(thisVol) ;        % Load Volume
tmp      = squeeze(Volume.vol) ;        % Squeeze useless dimensions
TR = Volume.pixdim(5) ;                 % Temporal Resolution

%% Slice Timing 

% Evenly space slice times between one TR
sliceSpace = 0:TR/length(uniqueTiming):TR;
sliceSpace = sliceSpace(1:end-1);

% Create vector with with each slice for every TR
sliceTimes = nan(length(sliceSpace),size(tmp,4));
for i = 1:size(tmp,4)
    sliceTimes(:,i) = sliceSpace + (TR * (i - 1));
end
sliceTimes = sort(sliceTimes(:));

% Obtain avg signal for each slice in V1
avgSlice = nan(length(sliceSpace),size(tmp,4));
for i = 1:length(uniqueTiming)
    % tmp(X,Y,Z,T) - note that X and Y will need to change, based on V1
    maskSlices = Mask.vol(:,:,uniqueTiming(i)==thisTiming) ;
    V1slice = find(abs(maskSlices) == 1);
    vSlices = tmp(:,:,uniqueTiming(i)==thisTiming,:) ;
    RvolSlices = reshape(vSlices,size(vSlices,1)*size(vSlices,2)*size(vSlices,3),size(vSlices,4)) ;
    [psc] = convert_to_psc(RvolSlices);     % tc is N x TR --> percent signal change
    avgSlice(i,:) = mean(psc(V1slice,:)) ;  % avg signal for each slice (9 slices)
end

% Rearrange avg signal per slice -- index every x slice value
iSlice = nan(1,size(avgSlice,2)*size(avgSlice,1));
ct = 0;
for k = 1:size(avgSlice,2)
    for i = 1:size(avgSlice,1)
        ct = ct + 1;
        iSlice(ct) = avgSlice(i,k) ;
    end
end

%% Get Event Times
eventTimes = cell(length(stimDirs),1);
% loop through stimulus dirs -- Get Event Start Times
for i = 1:length(stimDirs)
    currFile = fullfile(stimDir,stimDirs(i)) ;
    % Get Event Start Time - file name
    eventTimeData = dir(char(fullfile(currFile,eventTimesFile))) ;
    % Load File (3 colum format)
    eventTimesLoad = load(char(fullfile(currFile,eventTimeData.name))) ;
    % Save first column (start times)
    eventTimes{i} = eventTimesLoad(:,1) ;
end

%% Derive HRF in V1 voxels
% freeview
sampT         = round(TR/length(uniqueTiming));
runEventTimes = round(eventTimes{1}*1000) ;
inputTC       = iSlice'; % choose the input time-series

% Derive the RF in avg V1 slices
HRF = deriveHRF(inputTC,runEventTimes,sampT,HRFdur,numFreqs);
hrf(ii,:) = resample(HRF,1,TR) ;

fprintf('Just finished run #%d\n', ii)
end

%% Plot pretty things
SEHRF = std(hrf)./sqrt(size(hrf,1));
BOLDHRF = mean(hrf) ;
figure;
% We donwsample because TR & HRFdur are in ms resolution
errorbar((0:TR:HRFdur-TR)/1000,BOLDHRF,SEHRF,'LineWidth',2)
xlabel('Time/s'); ylabel('Signal'); set(gca,'FontSize',15);
title('HRF');
 % Compare mean HRF with individual runs
figure;
plot(BOLDHRF,'k','LineWidth',8);hold on;
plot(hrf');

% Yay! End of Script