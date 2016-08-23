function [packets] = makePackets(sessionDir,packetType,func)

%   Outputs a 'packets' cell array of structures, each cell containing:
%
%   Usage:
%   [packets] = makePackets(sessionDir,packetType,func)
%
%   Outputs:
%   stimulus.values         - M x N matrix modeling M stimulus events
%   stimulus.timebase       - 1 x N vector of stimulus times (msec)
%   stimulus.metaData       - structure with info about the stimulus
%
%   response.values         - 1 x TR vector of response values
%   response.timebase       - 1 x TR vector of response times (msec)
%   response.metaData       - structure with info about the response
%
%   HRF.values              - 1 x N vector of response values
%   HRF.timebase            - 1 x N vector of response times (msec)
%   HRF.metaData            - structure with info about the stimulus
%
%   metaData.projectName    - project name (e.g. 'MelanopsinMR');
%   metaData.subjectName    - subject name (e.g. 'HERO_asb1');
%   metaData.sessionDate    - session date (e.g. '041416');
%   metaData.stimulusFile   - fullfile(sessionDir,'MatFiles',matFiles{i});
%   metaData.responseFile   - fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
%
%   Written by Andrew S Bock Aug 2016

%% set defaults
if ~exist('packetType','var') || isempty(packetType)
    packetType      = 'V1';
end
if ~exist('func','var') || isempty(func)
    func            = 'wdrf.tf';
end
%% File / path defaults
anatFileName        = 'mh.areas.anat.vol.nii.gz';
boldOutName         = 'mh.areas.func.vol.nii.gz';
bbregName           = 'func_bbreg.dat';
% stimulus files
matDir              = fullfile(sessionDir,'MatFiles');
matFiles            = listdir(matDir,'files');
% response files
boldDirs            = find_bold(sessionDir);
% HRF defaults
attentionTaskNames  = {'MirrorsOffMaxLMS','MirrorsOffMaxMel'};
HRFdur              = 16000;
numFreqs            = HRFdur/1000;
%% Meta data
[subjectStr,sessionDate]        = fileparts(sessionDir);
[projectStr,subjectName]        = fileparts(subjectStr);
[~,projectName]                 = fileparts(projectStr);
for i = 1:length(boldDirs)
    metaData{i}.projectName     = projectName;
    metaData{i}.subjectName     = subjectName;
    metaData{i}.sessionDate     = sessionDate;
    metaData{i}.stimulusFile    = fullfile(matDir,matFiles{i});
    metaData{i}.responseFile    = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
end
%% Load in fMRI data
for i = 1:length(boldDirs)
    inFile          = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
    inData{i}       = load_nifti(inFile);
end
%% Stimulus
for i = 1:length(boldDirs)
    % Get bold data details
    TR                              = inData{i}.pixdim(5); % TR in msec
    numTRs                          = size(inData{i}.vol,4);
    runDur                          = TR * numTRs; % length of run (msec)
    stimulus{i}.timebase            = 0:runDur-1;
    zVect                           = zeros(1,runDur);
    % Load that .mat file produced by the stimulus computer
    stimulus{i}.metaData            = load(fullfile(matDir,matFiles{i}));
    for j = 1:size(stimulus{i}.metaData.params.responseStruct.events,2)
        % phase offset
        if ~isempty(stimulus{i}.metaData.params.thePhaseOffsetSec)
            phaseOffsetSec = stimulus{i}.metaData.params.thePhaseOffsetSec(...
                stimulus{i}.metaData.params.thePhaseIndices(j));
        else
            phaseOffsetSec = 0;
        end
        % start time
        startTime = stimulus{i}.metaData.params.responseStruct.events(j).tTrialStart - ...
            stimulus{i}.metaData.params.responseStruct.tBlockStart + phaseOffsetSec;
        % duration
        if isfield(stimulus{i}.metaData.params.responseStruct.events(1).describe.params,'stepTimeSec')
            durTime = stimulus{i}.metaData.params.responseStruct.events(j).describe.params.stepTimeSec + ...
                2*stimulus{i}.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
        else
            durTime = stimulus{i}.metaData.params.responseStruct.events(j).tTrialEnd - ...
                stimulus{i}.metaData.params.responseStruct.events(j).tTrialStart;
        end
        % stimulus window
        stimWindow                  = ceil((startTime*1000) : (startTime*1000 + ((durTime*1000)-1)));
        % Save the stimulus values
        thisStim                    = zVect;
        thisStim(stimWindow)        = 1;
        % cosine ramp onset
        if stimulus{1}.metaData.params.responseStruct.events(1).describe.params.cosineWindowIn
            winDur  = stimulus{1}.metaData.params.responseStruct.events(1).describe.params.cosineWindowDurationSecs;
            cosOn   = (cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2;
            thisStim(stimWindow(1:winDur*1000)) = cosOn;
        end
        % cosine ramp offset
        if stimulus{1}.metaData.params.responseStruct.events(1).describe.params.cosineWindowOut
            winDur  = stimulus{1}.metaData.params.responseStruct.events(1).describe.params.cosineWindowDurationSecs;
            cosOff   = fliplr((cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2);
            thisStim(stimWindow(end-((winDur*1000)-1):end)) = cosOff;
        end
        % trim stimulus
        thisStim                    = thisStim(1:runDur); % trim events past end of run (occurs for stimuli presented near the end of the run)
        % save stimulus values
        stimulus{i}.values(j,:)     = thisStim;
    end
end
%% ROI
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
%% Response
% Convert fMRI to percent signal change, then average
for i = 1:length(boldDirs)
    % Get bold data details
    TR                              = inData{i}.pixdim(5); % TR in msec
    numTRs                          = size(inData{i}.vol,4);
    runDur                          = TR * numTRs; % length of run (msec)
    thisVol                         = inData{i}.vol;
    % reshape if a 4D volume
    if length(size(thisVol))>2
        thisVol                     = reshape(thisVol,...
            size(thisVol,1)*size(thisVol,2)*size(thisVol,3),size(thisVol,4));
    end
    ROIvals                         = thisVol(ROI{i},:);
    pscVals                         = convert_to_psc(ROIvals);
    timeSeries{i}                   = nanmean(pscVals);
    response{i}.values              = timeSeries'; % could also use 'cleanData'
    response{i}.timebase            = 0:TR:runDur-1; % beginning of each TR (msec)
    response{i}.metaData.filename   = fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
    response{i}.metaData.packetType = packetType;
end
%% HRF
for i = 1:length(boldDirs)
    % Get the attention events
    ct = 0;
    attEvents = [];
    for j = 1:size(stimulus{i}.metaData.params.responseStruct.events,2)
        switch metaData{i}.projectName
            case 'MelanopsinMR'
                % Get the attention events
                if sum(strcmp(stimulus{i}.metaData.params.responseStruct.events(j).describe.direction,attentionTaskNames))
                    ct = ct + 1;
                    % Get the stimulus window
                    attEvents(ct) = stimulus{i}.metaData.params.responseStruct.events(j).tTrialStart - ...
                        stimulus{i}.metaData.params.responseStruct.tBlockStart + ...
                        stimulus{i}.metaData.params.thePhaseOffsetSec(stimulus{i}.metaData.params.thePhaseIndices(j));
                end
            case 'MOUNT_SINAI'
                if stimulus{i}.metaData.params.responseStruct.events(j).attentionTask.segmentFlag
                    ct = ct + 1;
                    attEvents(ct) = stimulus{i}.metaData.params.responseStruct.events(j).t(...
                        stimulus{i}.metaData.params.responseStruct.events(j).attentionTask.T == 1) - ...
                        stimulus{i}.metaData.params.responseStruct.tBlockStart;
                end
        end
    end
    eventTimes                  = round(attEvents*1000); % attention events (msec)
    sampT                       = inData{i}.pixdim(5); % TR in msec
    [allHRF(i,:)]               = deriveHRF(...
        timeSeries{i}',eventTimes,sampT,HRFdur,numFreqs);
end
HRF.values                      = mean(allHRF); % average HRFs across runs (individual HRFs are noisy)
HRF.timebase                    = 0:HRFdur-1;
HRF.metaData.packetType         = packetType;
%% Save outputs
for i = 1:length(boldDirs)
    packets{i}.stimulus     = stimulus{i};
    packets{i}.response     = response{i};
    packets{i}.HRF          = HRF;
    packets{i}.metaData     = metaData{i};
end