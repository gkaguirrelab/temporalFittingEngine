function [HRF] = subjHRF(params)

% Calculates the subject-specific HRF across all sessions within a project
% directory for a given subject
%
%   Usage:
%   [HRF] = subjHRF(params)
%
%   Inputs:
%   params.subjDir     - subject directory (e.g. '/data/jag/MELA/MelanopsinMR/HERO_asb1')
%   params.roiType     - type of ROI (default = 'V1')
%   params.func        - functional volume (default = 'wdrf.tf')
%   params.eccRange    - eccentricity range for ROI (only used if roiType = 'V1stim' or 'V2V3stim');
%
%   Outputs:
%   HRF.mean    - mean HRF (1 x N vector, msec resolution)
%   HRF.sem     - standard error of them mean (1 x vector, msec resolution)
%   HRF.numRuns - number of runs used to calculate the HRF
%
%   Written by Andrew S Bock Aug 2016

%% Set defaults
if ~isfield(params,'roiType');
    params.roiType                  = 'V1';
end
if ~isfield(params,'func');
    params.func                     = 'wdrf.tf';
end
if ~isfield(params,'eccRange');
    params.eccRange                 = [2.5 32]; % based on MaxMel data
end
% pull out params fields
subjDir                             = params.subjDir;
roiType                             = params.roiType;
func                                = params.func;
eccRange                            = params.eccRange ; % based on the MaxMel data
% HRF defaults
attentionTaskNames  = {'MirrorsOffMaxLMS','MirrorsOffMaxMel','MirrorsOffSplatterControl'};
HRFdur              = 16000;
numFreqs            = HRFdur/1000;
% anatomical to functional registration file
bbregName                           = 'func_bbreg.dat';
%% Get the session directories
sessDirs = listdir(subjDir,'dirs');

%% Get the mean HRF within each session directory
allHRF = cell(1,length(sessDirs));
for k = 1:length(sessDirs)
    % session directory paths
    sessionDir                      = fullfile(subjDir,sessDirs{k});
    boldDirs                        = find_bold(sessionDir);
    matDir                          = fullfile(sessionDir,'MatFiles');
    matFiles                        = listdir(matDir,'files');
    % Meta data
    [subjectStr,~]                  = fileparts(sessionDir);
    [projectStr,~]                  = fileparts(subjectStr);
    [~,projectName]                 = fileparts(projectStr);
    allHRF{k}                       = nan(length(boldDirs),HRFdur);
    for i = 1:length(boldDirs)
        % Load fMRI data
        inFile                      = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
        inData                      = load_nifti(inFile);
        % Get the ROI voxels
        switch roiType
            case {'V1' 'V2V3'}
                areaInName          = 'mh.areas.anat.vol.nii.gz';
                areaOutName         = 'mh.areas.func.vol.nii.gz';
                eccInName           = 'mh.ecc.anat.vol.nii.gz';
                eccOutName          = 'mh.ecc.func.vol.nii.gz';
            case 'LGN'
                areaInName          = 'mh.LGN.nii.gz';
                areaOutName         = 'mh.LGN.func.vol.nii.gz';
        end
        % Project areanat to func
        areaInFile                  = fullfile(sessionDir,'anat_templates',areaInName);
        bbregFile                   = fullfile(sessionDir,boldDirs{i},bbregName); % registration file
        areaOutFile                 = fullfile(sessionDir,boldDirs{i},areaOutName);
        system(['mri_vol2vol --mov ' fullfile(sessionDir,boldDirs{i},[func '.nii.gz']) ...
            ' --targ ' areaInFile ' --o ' areaOutFile ...
            ' --reg ' bbregFile ' --inv --nearest']);
        areaData                    = load_nifti(areaOutFile);
        % If 'V1stim' or 'V2V3stim', also project eccentricity file
        switch roiType
            case {'V1' 'V2V3'}
                eccInFile           = fullfile(sessionDir,'anat_templates',eccInName);
                bbregFile           = fullfile(sessionDir,boldDirs{i},bbregName); % registration file
                eccOutFile          = fullfile(sessionDir,boldDirs{i},eccOutName);
                system(['mri_vol2vol --mov ' fullfile(sessionDir,boldDirs{i},[func '.nii.gz']) ...
                    ' --targ ' eccInFile ' --o ' eccOutFile ...
                    ' --reg ' bbregFile ' --inv --nearest']);
                eccData             = load_nifti(eccOutFile);
        end
        switch roiType
            case 'LGN'
                ROI                 = find(abs(areaData.vol)==1);
            case 'V1'
                ROI                 = find(abs(areaData.vol)==1 & ...
                    eccData.vol>eccRange(1) & eccData.vol<eccRange(2));
            case 'V2V3'
                ROI                 = find((abs(areaData.vol)==2 | abs(areaData.vol)==3) & ...
                    eccData.vol>eccRange(1) & eccData.vol<eccRange(2));
        end
        % Get the response
        thisVol                     = inData.vol;
        % reshape if a 4D volume
        if length(size(thisVol))>2
            thisVol                 = reshape(thisVol,...
                size(thisVol,1)*size(thisVol,2)*size(thisVol,3),size(thisVol,4));
        end
        ROIvals                     = thisVol(ROI,:);
        pscVals                     = convert_to_psc(ROIvals);
        % Get the median percent signal change for the ROI
        timeSeries                  = nanmedian(pscVals);
        % Get stimulus data
        stimulus.metaData           = load(fullfile(matDir,matFiles{i}));
        % Get the attention events
        ct = 0;
        attEvents = [];
        for j = 1:size(stimulus.metaData.params.responseStruct.events,2)
            switch projectName
                case 'MelanopsinMR'
                    % Get the attention events
                    if sum(strcmp(stimulus.metaData.params.responseStruct.events(j).describe.direction,attentionTaskNames))
                        ct = ct + 1;
                        % Get the stimulus window
                        attEvents(ct) = stimulus.metaData.params.responseStruct.events(j).tTrialStart - ...
                            stimulus.metaData.params.responseStruct.tBlockStart + ...
                            stimulus.metaData.params.thePhaseOffsetSec(stimulus.metaData.params.thePhaseIndices(j));
                    end
                case 'MOUNT_SINAI'
                    if stimulus.metaData.params.responseStruct.events(j).attentionTask.segmentFlag
                        ct = ct + 1;
                        attEvents(ct) = stimulus.metaData.params.responseStruct.events(j).t(...
                            stimulus.metaData.params.responseStruct.events(j).attentionTask.T == 1) - ...
                            stimulus.metaData.params.responseStruct.tBlockStart;
                    end
            end
        end
        eventTimes                  = round(attEvents*1000); % attention events (msec)
        sampT                       = inData.pixdim(5); % TR in msec
        allHRF{k}(i,:)              = deriveHRF(timeSeries',eventTimes,sampT,HRFdur,numFreqs);
    end
end
%% Save HRF
ct = 0;
for i = 1:length(allHRF)
    for j = 1:size(allHRF{i},1)
        ct = ct + 1;
        tmpHRF(ct,:) = allHRF{i}(j,:);
    end
end
tmpmean                             = nanmean(tmpHRF);
HRF.mean                            = tmpmean - tmpmean(1);
HRF.sem                             = nanstd(tmpHRF) / sqrt(size(tmpHRF,1));
HRF.numRuns                         = size(tmpHRF,1);
HRF.metaData.roiType                = roiType;
HRF.metaData.eccRange               = eccRange;
HRF.metaData.func                   = func;
for k = 1:length(sessDirs)
    % session directory paths
    sessionDir                      = fullfile(subjDir,sessDirs{k});
    hrfDir                          = fullfile(sessionDir,'HRF');
    if ~exist(hrfDir,'dir')
        mkdir(hrfDir);
    end
    save(fullfile(hrfDir,[roiType '.mat']),'HRF');
end