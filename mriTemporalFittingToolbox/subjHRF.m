function [HRF] = subjHRF(subjDir,nameROI,func)

% Calculates the subject-specific HRF across all sessions within a project
% directory for a given subject
%
%   Written by Andrew S Bock Aug 2016

%% Set defaults
% ROI
if ~exist('nameROI','var') || isempty(nameROI)
    nameROI                         = 'V1';
end
% functional volume
if ~exist('func','var') || isempty(func)
    func                            = 'wdrf.tf';
end
% HRF defaults
attentionTaskNames  = {'MirrorsOffMaxLMS','MirrorsOffMaxMel','MirrorsOffSplatterControl'};
HRFdur              = 16000;
numFreqs            = HRFdur/1000;
% anatomical to functional registration file
bbregName                           = 'func_bbreg.dat';
%% Get the session directories
sessDirs = listdir(subjDir,'dirs');

%% Get the mean HRF within each session directory
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
    progBar = ProgressBar(length(boldDirs),'fooing...');
    for i = 1:length(boldDirs)
        % Load fMRI data
        inFile                      = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
        inData                      = load_nifti(inFile);
        % Get the ROI voxels
        switch nameROI
            case {'V1' 'V2V3'}
                anatFileName        = 'mh.areas.anat.vol.nii.gz';
                boldOutName         = 'mh.areas.func.vol.nii.gz';
            case 'LGN'
                anatFileName        = 'mh.LGN.nii.gz';
                boldOutName         = 'mh.LGN.func.vol.nii.gz';
        end
        anatFile                    = fullfile(sessionDir,'anat_templates',anatFileName);
        bbregFile                   = fullfile(sessionDir,boldDirs{i},bbregName); % registration file
        outFile                     = fullfile(sessionDir,boldDirs{i},boldOutName);
        system(['mri_vol2vol --mov ' fullfile(sessionDir,boldDirs{i},[func '.nii.gz']) ...
            ' --targ ' anatFile ' --o ' outFile ...
            ' --reg ' bbregFile ' --inv --nearest']);
        areaData                    = load_nifti(outFile);
        switch nameROI
            case {'V1' 'LGN'}
                ROI                 = find(abs(areaData.vol)==1);
            case 'V2V3'
                ROI                 = find(abs(areaData.vol)==2 | abs(areaData.vol)==3);
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
        timeSeries                  = nanmean(pscVals);    
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
        [allHRF(k,i,:)]             = deriveHRF(timeSeries',eventTimes,sampT,HRFdur,numFreqs);
        progBar(i);
    end
end