function packet = makePacket(params)

%   Outputs a 'packet' structure with stimulus, response, metaData, and
%   (optionally) HRF information
%
%   Usage:
%   packet = makePacket(params)
%
%   Input params structure:
%   params.sessionDir
%   params.runNum
%   params.responseFile
%   params.stimulusFile
%   params.hrfFile
%   params.timeSeries
%   params.TR
%   params.roiType
%   params.packetType
%
%   Output fields in packets:
%
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
%
%   HRF.values              - 1 x N vector of response values
%   HRF.timebase            - 1 x N vector of response times (msec)
%   HRF.metaData            - structure with info about the HRF
%
%   Written by Andrew S Bock Aug 2016

%% Set defaults
switch params.packetType
    case 'bold'
        runNames                            = find_bold(params.sessionDir);
    case 'pupil'
        runNames                            = listdir(fullfile(params.sessionDir, 'EyeTrackingFiles/*.mat'), 'files');
        params.LiveTrackSamplingRate        = 60; % Hz
        params.ResamplingFineFreq           = 1000; % 1 msec
        params.BlinkWindowSample            = -50:50; % Samples surrounding the blink event
        params.TRDurSecs                    = 0.8;
end
if isempty(runNames)
    error(['No runs found in ' sessionDir]);
end
% stimulus files
%% Metadata
[subjectStr,sessionDate]                    = fileparts(params.sessionDir);
[projectStr,subjectName]                    = fileparts(subjectStr);
[~,projectName]                             = fileparts(projectStr);
metaData.projectName                        = projectName;
metaData.subjectName                        = subjectName;
metaData.sessionDate                        = sessionDate;
metaData.stimulusFile                       = params.stimulusFile;
metaData.responseFile                       = params.responseFile;
%% Stimulus
% Load that .mat file produced by the stimulus computer
stimulus.metaData                           = load(params.stimulusFile);
% Get run duration
runDur                                      = sum(stimulus.metaData.params.trialDuration)*1000; % length of run (msec)
% Set the timebase
stimulus.timebase                           = 0:runDur-1;
zVect                                       = zeros(1,runDur);
for j = 1:size(stimulus.metaData.params.responseStruct.events,2)
    % phase offset
    if ~isempty(stimulus.metaData.params.thePhaseOffsetSec)
        phaseOffsetSec = stimulus.metaData.params.thePhaseOffsetSec(...
            stimulus.metaData.params.thePhaseIndices(j));
    else
        phaseOffsetSec = 0;
    end
    % start time
    startTime = stimulus.metaData.params.responseStruct.events(j).tTrialStart - ...
        stimulus.metaData.params.responseStruct.tBlockStart + phaseOffsetSec;
    % duration
    if isfield(stimulus.metaData.params.responseStruct.events(1).describe.params,'stepTimeSec')
        durTime = stimulus.metaData.params.responseStruct.events(j).describe.params.stepTimeSec + ...
            2*stimulus.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
    else
        durTime = stimulus.metaData.params.responseStruct.events(j).tTrialEnd - ...
            stimulus.metaData.params.responseStruct.events(j).tTrialStart;
    end
    % stimulus window
    stimWindow                              = ceil((startTime*1000) : (startTime*1000 + ((durTime*1000)-1)));
    % Save the stimulus values
    thisStim                                = zVect;
    thisStim(stimWindow)                    = 1;
    % cosine ramp onset
    if stimulus.metaData.params.responseStruct.events(j).describe.params.cosineWindowIn
        winDur  = stimulus.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
        cosOn   = (cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2;
        thisStim(stimWindow(1:winDur*1000)) = cosOn;
    end
    % cosine ramp offset
    if stimulus.metaData.params.responseStruct.events(j).describe.params.cosineWindowOut
        winDur  = stimulus.metaData.params.responseStruct.events(j).describe.params.cosineWindowDurationSecs;
        cosOff   = fliplr((cos(pi+linspace(0,1,winDur*1000)*pi)+1)/2);
        thisStim(stimWindow(end-((winDur*1000)-1):end)) = cosOff;
    end
    % trim stimulus
    thisStim                                = thisStim(1:runDur); % trim events past end of run (occurs for stimuli presented near the end of the run)
    % save stimulus values
    stimulus.values(j,:)                    = thisStim;
end
%% Response
switch params.packetType
    case 'bold'
        % Get bold data details
        numTRs                              = size(params.timeSeries,2);
        respDur                             = params.TR * numTRs; % length of run (msec)
        response.values                     = timeSeries'; % could also use 'cleanData'
        response.timebase                   = 0:params.TR:respDur-1; % beginning of each TR (msec)
        response.metaData.filename          = params.responseFile;
        response.metaData.roiType           = params.roiType;
        % HRF (if applicable)
        switch roiType
            case {'LGN' 'V1' 'V2V3'}
                tmp                         = load(params.hrfFile);
                HRF.values                  = tmp.mean;
                HRF.timebase                = 0:length(HRF.values)-1;
                HRF.metaData                = tmp.metaData;
        end
    case 'pupil'
        response.timebase                   = stimulus.timebase;
        params.TimeVectorFine               = response.timebase;
        switch sessionDate
            case {'053116' '060116' '060216'}
                params.acquisitionFreq      = 30;
            otherwise
                params.acquisitionFreq      = 60;
        end
        params.NTRsExpected                 = runDur/(params.TRDurSecs*1000);
        response.values                     = loadPupilDataForPackets(fullfile(params.sessionDir, 'EyeTrackingFiles', runNames{params.runNum}), stimulus, metaData, params);
end
%% Save the packets
packet.stimulus                     = stimulus;
packet.response                     = response;
switch packetType
    case 'bold'
        packet.HRF                  = HRF;
end
packet.metaData                     = metaData;