function [packets] = makePackets(sessionDir,packetType,func)

%   Outputs a 'packets' cell array of structures, each cell containing:
%
%   stimulus.values     - M x N matrix modeling M stimulus events
%   stimulus.timebase   - 1 x N vector of stimulus times (msec)
%   stimulus.file       - 1 x M cell array with the full path to the stimulus text file
%
%   response.values     - 1 x TR vector of response values
%   response.timebase   - 1 x TR vector of response times (msec)
%
%   HRF.values          - 1 x N vector of response values
%   HRF.timebase        - 1 x N vector of response times (msec, should start at 0)
%
%   Usage:
%   [packets] = makePackets(inFile,packetType,anatFile,stimFile)
%
%   Written by Andrew S Bock Aug 2016

%% set defaults
if ~exist('packetType','var') || isempty(packetType)
    packetType      = 'V1';
end
if ~exist('func','var') || isempty(func)
    func            = 'wdrf.tf';
end
% File / path defaults
anatFileName        = 'mh.areas.anat.vol.nii.gz';
boldOutName         = 'mh.areas.func.vol.nii.gz';
bbregName           = 'func_bbreg.dat';
boldDirs            = find_bold(sessionDir);
stimDirs            = listdir(fullfile(sessionDir,'Stimuli'),'dirs');
% HRF defaults
HRFdur              = 16000;
numFreqs            = HRFdur/1000;
%% Load in fMRI data
for i = 1:length(boldDirs)
    inFile          = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
    inData{i}       = load_nifti(inFile);
end
%% Stimulus
for i = 1:length(boldDirs)
    % Get bold data details
    TR      = inData{i}.pixdim(5); % TR in msec
    numTRs  = size(inData{i}.vol,4);
    runDur  = TR * numTRs; % length of run (msec)
    zVect   = zeros(1,runDur);
    ct = 0;
    % stimulus files
    stimFiles = listdir(fullfile(sessionDir,'Stimuli',stimDirs{i},'*_valid.txt'),'files');
    for j = 1:length(stimFiles)
        % different stimulus types are specified by different stimulus files
        stimData = load(fullfile(sessionDir,'Stimuli',stimDirs{i},stimFiles{j}));
        % stimulus events for each stimulus type
        for k = 1:size(stimData,1)
            ct = ct + 1;
            tmpWindow                   = ...
                (stimData(k,1)*1000) : ( (stimData(k,1)*1000) + (stimData(k,2)*1000)-1 );
            tmpWindow                   = ceil(tmpWindow); % for event timing greater than msec precision, typically starting ~0 seconds
            thisVol                     = zVect;
            thisVol(tmpWindow)          = 1;
            thisVol                     = thisVol(1:runDur); % trim events past end of run (occurs for stimuli presented near the end of the run)
            stimulus{i}.values(ct,:)    = thisVol;
            stimulus{i}.timebase(ct,:)  = 0:TR:(runDur - TR);
            stimulus{i}.file{ct}        = fullfile(sessionDir,'Stimuli',stimDirs{i},stimFiles{j});
        end
    end
end
%% Response
switch packetType
    case 'V1'
        anatFile        = fullfile(sessionDir,'anat_templates',anatFileName);
        for i = 1:length(boldDirs)
            % project anatomical file to functional space for each run
            bbregFile   = fullfile(sessionDir,boldDirs{i},bbregName); % registration file
            outFile     = fullfile(sessionDir,boldDirs{i},boldOutName);
            system(['mri_vol2vol --mov ' fullfile(sessionDir,boldDirs{i},[func '.nii.gz']) ...
                ' --targ ' anatFile ' --o ' outFile ...
                ' --reg ' bbregFile ' --inv --nearest']);
            areaData    = load_nifti(outFile);
            ROI{i}      = find(abs(areaData.vol)==1);
        end
end
% Convert fMRI to percent signal change, then average
for i = 1:length(boldDirs)
    % Get bold data details
    TR                  = inData{i}.pixdim(5); % TR in msec
    numTRs              = size(inData{i}.vol,4);
    runDur              = TR * numTRs; % length of run (msec)
    thisVol             = inData{i}.vol;
    % reshape if a 4D volume
    if length(size(thisVol))>2
        thisVol         = reshape(thisVol,size(thisVol,1)*size(thisVol,2)*size(thisVol,3),size(thisVol,4));
    end
    ROIvals             = thisVol(ROI{i},:);
    pscVals{i}          = convert_to_psc(ROIvals);
    meanVals{i}         = nanmean(pscVals{i});
end
%% HRF
for i = 1:length(boldDirs)
    % Get bold data details
    timeSeries              = meanVals{i}'; % TR x N
    % Get the attention events
    attFile                 = listdir(fullfile(sessionDir,'Stimuli',stimDirs{i},...
        '*_attentionTask.txt'),'files');
    attEvents               = load(fullfile(sessionDir,'Stimuli',stimDirs{i},...
        attFile{1}));
    eventTimes              = round(attEvents(:,1)*1000); % attention events (msec)
    sampT                   = inData{i}.pixdim(5); % TR in msec
    [allHRF(i,:),cleanData]   = deriveHRF(...
        timeSeries,eventTimes,sampT,HRFdur,numFreqs);
    response{i}.values      = cleanData';  
    response{i}.timebase    = 0:TR:(runDur - TR);
end
HRF.values                  = mean(allHRF); % average HRFs across runs (individual HRFs are noisy)
HRF.timebase                = 0:HRFdur-1;
%% Save outputs
for i = 1:length(boldDirs)
    packets{i}.stimulus     = stimulus{i};
    packets{i}.response     = response{i};
    packets{i}.HRF          = HRF; 
end