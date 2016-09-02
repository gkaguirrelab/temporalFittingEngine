function packets = makePackets(sessionDir,packetType,roiType,func,saveFlag)

%   Outputs a 'packets' cell array of structures, each cell containing:
%
%   Usage:
%   packets = makePackets(sessionDir,packetType,roiType,func)
%
%   Output fields in packets:
%   stimulus.values         - M x N matrix modeling M stimulus events
%   stimulus.timebase       - 1 x N vector of stimulus times (msec)
%   stimulus.metaData       - structure with info about the stimulus
%
%   response.values         - 1 x TR vector of response values
%   response.timebase       - 1 x TR vector of response times (msec)
%   response.metaData       - structure with info about the response
%
%   metaData.projectName    - project name (e.g. 'MelanopsinMR');
%   metaData.subjectName    - subject name (e.g. 'HERO_asb1');
%   metaData.sessionDate    - session date (e.g. '041416');
%   metaData.stimulusFile   - fullfile(sessionDir,'MatFiles',matFiles{i});
%   metaData.responseFile   - fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
%
%   If packetType == 'bold', also outputs:
%   HRF.values              - 1 x N vector of response values
%   HRF.timebase            - 1 x N vector of response times (msec)
%   HRF.metaData            - structure with info about the HRF
%
%   Written by Andrew S Bock Aug 2016

%% Set defaults
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = true;
end
switch packetType
    case 'bold'
        % set defaults
        if ~exist('roiType','var') || isempty(roiType)
            roiType                         = 'V1';
        end
        if ~exist('func','var') || isempty(func)
            func                            = 'wdrf.tf';
        end
        % File / path defaults
        switch roiType
            case {'V1' 'V2V3'}
                anatFileName                = 'mh.areas.anat.vol.nii.gz';
                boldOutName                 = 'mh.areas.func.vol.nii.gz';
            case 'LGN'
                anatFileName                = 'mh.LGN.nii.gz';
                boldOutName                 = 'mh.LGN.func.vol.nii.gz';
        end
        % anatomical to functional registration file
        bbregName                           = 'func_bbreg.dat';
        % HRF defaults
        hrfDir                              = fullfile(sessionDir,'HRF');
        % response files
        runNames                            = find_bold(sessionDir);
    case 'pupil'
        runNames = listdir(fullfile(sessionDir, 'EyeTrackingFiles/*.mat'), 'files');
        params.LiveTrackSamplingRate        = 60; % Hz
        params.ResamplingFineFreq           = 1000; % 1 msec
        params.BlinkWindowSample            = -50:50; % Samples surrounding the blink event
        params.TRDurSecs                    = 0.8;
end
if isempty(runNames)
   error(['No runs found in ' sessionDir]); 
end
% stimulus files
matDir                                      = fullfile(sessionDir,'MatFiles');
matFiles                                    = listdir(matDir,'files');
if saveFlag
    % save directory
    saveDir                                     = fullfile(sessionDir,'Packets');
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
end
%% Metadata
[subjectStr,sessionDate]                    = fileparts(sessionDir);
[projectStr,subjectName]                    = fileparts(subjectStr);
[~,projectName]                             = fileparts(projectStr);
for i = 1:length(runNames)
    metaData{i}.projectName                 = projectName;
    metaData{i}.subjectName                 = subjectName;
    metaData{i}.sessionDate                 = sessionDate;
    metaData{i}.stimulusFile                = fullfile(matDir,matFiles{i});
    switch packetType
        case 'bold'
            metaData{i}.responseFile        = fullfile(sessionDir,runNames{i},[func '.nii.gz']);
        case 'pupil'
            metaData{i}.responseFile        = fullfile(sessionDir,runNames{i});
    end
end
%% Stimulus
for i = 1:length(runNames)
    % Load that .mat file produced by the stimulus computer
    stimulus{i}.metaData                    = load(fullfile(matDir,matFiles{i}));
    % Get run duration
    runDur                                  = sum(stimulus{i}.metaData.params.trialDuration)*1000; % length of run (msec)
    % Set the timebase
    stimulus{i}.timebase                    = 0:runDur-1;
    zVect                                   = zeros(1,runDur);
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
        stimWindow                          = ceil((startTime*1000) : (startTime*1000 + ((durTime*1000)-1)));
        % Save the stimulus values
        thisStim                            = zVect;
        thisStim(stimWindow)                = 1;
        % cosine ramp onset
        if stimulus{i}.metaData.params.responseStruct.events(j).describe.params.cosineWindowIn
            winDur  = stimulus{i}.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
            cosOn   = (cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2;
            thisStim(stimWindow(1:winDur*1000)) = cosOn;
        end
        % cosine ramp offset
        if stimulus{i}.metaData.params.responseStruct.events(j).describe.params.cosineWindowOut
            winDur  = stimulus{i}.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
            cosOff   = fliplr((cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2);
            thisStim(stimWindow(end-((winDur*1000)-1):end)) = cosOff;
        end
        % trim stimulus
        thisStim                            = thisStim(1:runDur); % trim events past end of run (occurs for stimuli presented near the end of the run)
        % save stimulus values
        stimulus{i}.values(j,:)             = thisStim;
    end
end
%% Response
switch packetType
    case 'bold'
        % Load in fMRI data
        for i = 1:length(runNames)
            inFile                          = fullfile(sessionDir,runNames{i},[func '.nii.gz']);
            inData{i}                       = load_nifti(inFile);
        end
        % ROI
        anatFile            = fullfile(sessionDir,'anat_templates',anatFileName);
        for i = 1:length(runNames)
            % project anatomical file to functional space for each run
            bbregFile       = fullfile(sessionDir,runNames{i},bbregName); % registration file
            outFile         = fullfile(sessionDir,runNames{i},boldOutName);
            system(['mri_vol2vol --mov ' fullfile(sessionDir,runNames{i},[func '.nii.gz']) ...
                ' --targ ' anatFile ' --o ' outFile ...
                ' --reg ' bbregFile ' --inv --nearest']);
            areaData        = load_nifti(outFile);
            switch roiType
                case {'V1' 'LGN'}
                    ROI{i}  = find(abs(areaData.vol)==1);
                case 'V2V3'
                    ROI{i}  = find(abs(areaData.vol)==2 | abs(areaData.vol)==3);
            end
        end
        % Convert fMRI to percent signal change, then average
        for i = 1:length(runNames)
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
            response{i}.metaData.filename   = fullfile(sessionDir,runNames{i},[func '.nii.gz']);
            response{i}.metaData.roiType = roiType;
        end
        % HRF (if applicable)
        switch roiType
            case {'LGN' 'V1' 'V2V3'}
                tmp                         = load(fullfile(hrfDir,[roiType '.mat']));
                HRF.values                  = tmp.mean;
                HRF.timebase                = 0:length(HRF.values)-1;
                HRF.metaData                = tmp.metaData;
        end
    case 'pupil'
        for i = 1:length(runNames)
            response{i}.timebase            = stimulus{i}.timebase;
            params.TimeVectorFine           = response{i}.timebase;
            switch sessionDate
                case {'053116' '060116' '060216'}
                    params.acquisitionFreq  = 30;
                otherwise
                    params.acquisitionFreq  = 60;
            end
            params.NTRsExpected             = runDur/(params.TRDurSecs*1000);
            response{i}.values = loadPupilDataForPackets(fullfile(sessionDir, 'EyeTrackingFiles', runNames{i}), stimulus{i}, metaData{i}, params);
        end
end
%% Save the packets
for i = 1:length(runNames)
    packets{i}.stimulus                     = stimulus{i};
    packets{i}.response                     = response{i};
    switch packetType
        case 'bold'
            packets{i}.HRF                  = HRF;
    end
    packets{i}.metaData                     = metaData{i};
end
if saveFlag
    save(fullfile(saveDir,[roiType '.mat']),'packets','-v7.3');
end