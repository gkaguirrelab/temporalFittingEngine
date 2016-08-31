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
                    acquisitionFreq = 60;
                otherwise
                    acquisitionFreq = 30;
            end
            
            %
            %%% EYE TRACKING DATA %%%
            % Load the eye tracking data
            Data_LiveTrack = load(fullfile(sessionDir, 'EyeTrackingFiles', runDirs{i}));
            
            % Find cases in which the TTL pulse signal was split over the two
            % samples, and remove the second sample.
            Data_LiveTrack_TTLPulses_raw = [Data_LiveTrack.params.Report.Digital_IO1];
            tmpIdx = strfind(Data_LiveTrack_TTLPulses_raw, [1 1]);
            Data_LiveTrack_TTLPulses_raw(tmpIdx) = 0;
            
            Data_LiveTrack_TTLPulses = [];
            Data_LiveTrack_PupilDiameter = [];
            Data_LiveTrack_IsTracked = [];
            % We reconstruct the data set collected at 30/60 Hz.
            for rr = 1:length(Data_LiveTrack.params.Report)
                % Depending on how we set up the acquisiton frequency, we have to
                % do different things to extract the data
                switch acquisitionFreq
                    case 60
                        Data_LiveTrack_TTLPulses = [Data_LiveTrack_TTLPulses Data_LiveTrack_TTLPulses_raw(rr) 0];
                        
                        % We use the pupil width as the index of pupil diameter
                        Data_LiveTrack_PupilDiameter = [Data_LiveTrack_PupilDiameter Data_LiveTrack.params.Report(rr).PupilWidth_Ch01 ...
                            Data_LiveTrack.params.Report(rr).PupilWidth_Ch02];
                        
                        % Special case
                        if strcmp(dateID, '060616') && strcmp(subjectID, 'HERO_gka1');
                            Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked ...
                                Data_LiveTrack.params.Report(rr).S2Tracked];
                        else
                            Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked_Ch01 ...
                                Data_LiveTrack.params.Report(rr).PupilTracked_Ch02];
                        end
                    case 30
                        Data_LiveTrack_TTLPulses = [Data_LiveTrack_TTLPulses Data_LiveTrack_TTLPulses_raw(rr)];
                        Data_LiveTrack_PupilDiameter = [Data_LiveTrack_PupilDiameter Data_LiveTrack.params.Report(rr).LeftPupilWidth];
                        Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked];
                        
                end
            end
            
            % Now, we reconstruct the time vector of the data.
            TTLPulseIndices = find(Data_LiveTrack_TTLPulses); FirstTTLPulse = TTLPulseIndices(1);
            TimeVectorLinear = zeros(1, size(Data_LiveTrack_TTLPulses, 2));
            TimeVectorLinear(TTLPulseIndices) = (1:params.NTRsExpected)-1;
            
            % Replace zeros with NaN
            tmpX = 1:length(TimeVectorLinear);
            TimeVectorLinear(TimeVectorLinear == 0) = NaN;
            TimeVectorLinear(isnan(TimeVectorLinear)) = interp1(tmpX(~isnan(TimeVectorLinear)), ...
                TimeVectorLinear(~isnan(TimeVectorLinear)), tmpX(isnan(TimeVectorLinear)), 'linear', 'extrap');
            
            % Resample the timing to 1 msecs sampling
            Data_LiveTrack_PupilDiameter_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
                Data_LiveTrack_PupilDiameter, params.TimeVectorFine);
            Data_LiveTrack_IsTracked_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
                Data_LiveTrack_IsTracked, params.TimeVectorFine, 'nearest'); % Use NN interpolation for the binary tracking state

            % Extract the stimulus timing
            keyPressWhich = [];
            keyPressWhen = [];
            for rr = 1:length(stimulus{i}.metaData.params.responseStruct.events)
                keyPressWhich = [keyPressWhich stimulus{i}.metaData.params.responseStruct.events(rr).buffer.keyCode];
                keyPressWhen = [keyPressWhen stimulus{i}.metaData.params.responseStruct.events(rr).buffer.when];
            end
            
            % Extract only the ts
            keyPressWhen = keyPressWhen(keyPressWhich == 18);
            
            % Tack the first t also in this vector
            stimulus{i}.metaData_TTL = [stimulus{i}.metaData.params.responseStruct.tBlockStart keyPressWhen];
            
            % Subtract the absolute time of the first t
            stimulus{i}.metaData_TTL_t0 = stimulus{i}.metaData_TTL(1);
            stimulus{i}.metaData_TTL = stimulus{i}.metaData_TTL-stimulus{i}.metaData_TTL_t0;
            
            % Check that we have as many TRs as we expect
            fprintf('> Expecting <strong>%g</strong> TRs - Found <strong>%g</strong> (LiveTrack) and <strong>%g</strong> (OneLight record).\n', ...
                params.NTRsExpected, sum(Data_LiveTrack_TTLPulses), length(stimulus{i}.metaData_TTL));
            if (params.NTRsExpected == sum(Data_LiveTrack_TTLPulses)) || (params.NTRsExpected == length(stimulus{i}.metaData_TTL))
                fprintf('\t>> Expected number of TRs matches actual number.\n');
            else
                error('\t>> Mismatch between expected and actual number of TRs received.');
            end
            
            % We now have four variables of interest
            %   Data_LiveTrack_IsTracked_FineMasterTime <- Binary array indicating tracking state
            %   Data_LiveTrack_PupilDiameter_FineMasterTime <- Pupil diameter
            %   params.TimeVectorFine <- Time vector in TR time
            
            % Remove blinks from the pupil data
            Data_LiveTrack_BlinkIdx = [];
            for rr = 1:length(params.BlinkWindowSample)
                Data_LiveTrack_BlinkIdx = [Data_LiveTrack_BlinkIdx find(~Data_LiveTrack_IsTracked_FineMasterTime)+params.BlinkWindowSample(rr)];
            end
            % Remove any blinks from before the first sample
            Data_LiveTrack_BlinkIdx(Data_LiveTrack_BlinkIdx < 1) = []; %% ALSO CUT THE BLINKING OFF AT THE END
            
            % Remove any blinks after the last sample
            Data_LiveTrack_BlinkIdx(Data_LiveTrack_BlinkIdx > length(Data_LiveTrack_IsTracked_FineMasterTime)) = []; %% ALSO CUT THE BLINKING OFF AT THE END
            Data_LiveTrack_BlinkIdx = unique(Data_LiveTrack_BlinkIdx);
            Data_LiveTrack_PupilDiameter_FineMasterTime(Data_LiveTrack_BlinkIdx) = NaN;
            
            % Interpolate the elements
            Data_LiveTrack_PupilDiameter_FineMasterTime(isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)) = interp1(params.TimeVectorFine(~isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)), Data_LiveTrack_PupilDiameter_FineMasterTime(~isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)), params.TimeVectorFine(isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)));
            % Get the DC
            Data_LiveTrack_PupilDiameter_FineMasterTime_DC = nanmean(Data_LiveTrack_PupilDiameter_FineMasterTime);
            
            Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered = (Data_LiveTrack_PupilDiameter_FineMasterTime-nanmean(Data_LiveTrack_PupilDiameter_FineMasterTime))./(nanmean(Data_LiveTrack_PupilDiameter_FineMasterTime));
            
            
            % Low-pass filter the pupil data
            % Set up filter properties
            NFreqsToFilter = 8; % Number of low frequencies to remove
            for ii = 1:NFreqsToFilter
                X(2*ii-1,:) = sin(linspace(0, 2*pi*ii, size(stimulus{i}.timebase, 2)));
                X(2*ii,:) = cos(linspace(0, 2*pi*ii, size(stimulus{i}.timebase, 2)));
            end
            
            % Filter it
            [b, bint, r] = regress(Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered',X');
                subplot(3,1,1)
                plot(params.TimeVectorFine, Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered)
                subplot(3,1,2)
                plot(params.TimeVectorFine, r)
                subplot(3,1,3)
                plot(params.TimeVectorFine, Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered-r')
            
            % Create the filtered version
            Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered_f = r';
            
            % Add the DC back in
            Data_LiveTrack_PupilDiameter_FineMasterTime_f = Data_LiveTrack_PupilDiameter_FineMasterTime_DC + ...
                Data_LiveTrack_PupilDiameter_FineMasterTime_DC*Data_LiveTrack_PupilDiameter_FineMasterTime_MeanCentered_f;
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
        %% Do a bunch of pupil-related things
        %
        
        %% Save outputs
        for i = 1:length(runDirs)
            packets{i}.stimulus     = stimulus{i};
            packets{i}.response     = response{i};
            packets{i}.metaData     = metaData{i};
        end
        save(fullfile(saveDir,'.mat'),'packets','-v7.3');
end