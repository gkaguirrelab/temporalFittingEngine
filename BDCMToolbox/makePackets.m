function packets = makePackets(sessionDir,packetType,roiType,func)

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
%   Written by Andrew S Bock Aug 2016

switch packetType
    case 'bold'
        %% set defaults
        if ~exist('roiType','var') || isempty(roiType)
            roiType      = 'V1';
        end
        if ~exist('func','var') || isempty(func)
            func            = 'wdrf.tf';
        end
        %% File / path defaults
        switch roiType
            case {'V1' 'V2V3'}
                anatFileName        = 'mh.areas.anat.vol.nii.gz';
                boldOutName         = 'mh.areas.func.vol.nii.gz';
            case 'LGN'
                anatFileName        = 'mh.LGN.nii.gz';
                boldOutName         = 'mh.LGN.func.vol.nii.gz';
        end
        % anatomical to functional registration file
        bbregName           = 'func_bbreg.dat';
        % HRF defaults
        hrfDir              = fullfile(sessionDir,'HRF');
        % response files
        runDirs             = find_bold(sessionDir);
    case 'pupil'
        runDirs = listdir(fullfile(sessionDir, 'EyeTrackingFiles/*.mat'), 'files');
end

% stimulus files
matDir              = fullfile(sessionDir,'MatFiles');
matFiles            = listdir(matDir,'files');

% save directory
saveDir             = fullfile(sessionDir,'Packets');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Meta data
[subjectStr,sessionDate]        = fileparts(sessionDir);
[projectStr,subjectName]        = fileparts(subjectStr);
[~,projectName]                 = fileparts(projectStr);
for i = 1:length(runDirs)
    metaData{i}.projectName     = projectName;
    metaData{i}.subjectName     = subjectName;
    metaData{i}.sessionDate     = sessionDate;
    metaData{i}.stimulusFile    = fullfile(matDir,matFiles{i});
    switch packetType
        case 'bold'
            metaData{i}.responseFile    = fullfile(sessionDir,runDirs{i},[func '.nii.gz']);
        case 'pupil'
            metaData{i}.responseFile = fullfile(sessionDir,runDirs{i});
    end
end


%% Response data loading
switch packetType
    case 'bold'
        %% Load in fMRI data
        for i = 1:length(runDirs)
            inFile          = fullfile(sessionDir,runDirs{i},[func '.nii.gz']);
            inData{i}       = load_nifti(inFile);
        end
        
        %% ROI
        anatFile            = fullfile(sessionDir,'anat_templates',anatFileName);
        for i = 1:length(runDirs)
            % project anatomical file to functional space for each run
            bbregFile       = fullfile(sessionDir,runDirs{i},bbregName); % registration file
            outFile         = fullfile(sessionDir,runDirs{i},boldOutName);
            system(['mri_vol2vol --mov ' fullfile(sessionDir,runDirs{i},[func '.nii.gz']) ...
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
        
        %% Response
        % Convert fMRI to percent signal change, then average
        for i = 1:length(runDirs)
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
            response{i}.metaData.filename   = fullfile(sessionDir,runDirs{i},[func '.nii.gz']);
            response{i}.metaData.roiType = roiType;
        end
        
        %% HRF (if applicable)
        switch roiType
            case {'V1' 'LGN' 'V2V3'}
                HRF                         = load(fullfile(hrfDir,[roiType '.mat']));
        end
end


%% Stimulus
for i = 1:length(runDirs)
    % Load that .mat file produced by the stimulus computer
    stimulus{i}.metaData            = load(fullfile(matDir,matFiles{i}));
    
    % Get bold data details
    switch packetType
        case 'bold'
            TR                              = inData{i}.pixdim(5); % TR in msec
            numTRs                          = size(inData{i}.vol,4);
            runDur                          = TR * numTRs;
        case 'pupil'
            % Extract the run duration from the stimulus files since we do
            % not have any other information available for the pupil data.
            runDur                          = sum(stimulus{1}.metaData.params.trialDuration)*1000; % length of run (msec)
    end
    
    stimulus{i}.timebase            = 0:runDur-1;
    zVect                           = zeros(1,runDur);
    
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
        thisStim                    = thisStim(1:runDur); % trim events past end of run (occurs for stimuli presented near the end of the run)
        % save stimulus values
        stimulus{i}.values(j,:)     = thisStim;
    end
end

%% Response packet loading pupil
switch packetType
    case 'pupil'
        params.LiveTrackSamplingRate = 60; % Hz
        params.ResamplingFineFreq = 1000; % 1 msec
        params.BlinkWindowSample = -50:50; % Samples surrounding the blink event
        params.NTRsExpected = runDur/800;
        params.TRDurSecs = 0.8;
        
        for i = 1:length(runDirs)
            response{i}.timebase = stimulus{i}.timebase;
            params.TimeVectorFine = response{i}.timebase;
            
            switch sessionDate
                case {'053116' '060116' '060216'}
                    params.acquisitionFreq = 60;
                otherwise
                    params.acquisitionFreq = 30;
            end
            response{i}.values = loadPupilDataForPackets(fullfile(sessionDir, 'EyeTrackingFiles', runDirs{i}), stimulus, params);
        end
end

%% Saving data
switch packetType
    case 'bold'
        %% Save outputs
        for i = 1:length(runDirs)
            packets{i}.stimulus     = stimulus{i};
            packets{i}.response     = response{i};
            packets{i}.HRF          = HRF;
            packets{i}.metaData     = metaData{i};
        end
        save(fullfile(saveDir,[roiType '.mat']),'packets','-v7.3');
        
    case 'pupil'
        %% Save outputs
        for i = 1:length(runDirs)
            packets{i}.stimulus     = stimulus{i};
            packets{i}.response     = response{i};
            packets{i}.metaData     = metaData{i};
        end
        save(fullfile(saveDir,'.mat'),'packets','-v7.3');
end