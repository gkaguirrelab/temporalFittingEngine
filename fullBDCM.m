% This Script creates Temporal Transfer Functions 
% using a linear model and FIR with cosine windows.

%% Make sure toolboxes we need for this are on the path
AddToMatlabPathDynamically('BCDMToolbox');

%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_gka1' ; 
% *** Subject Pool ***
%     'HERO_asb1'--NOTE THAT ASB1 HAS 2 EXTRA LIGHTFLUX B AND 2 EXTRA S A
%                  RUNS--WE USE THE FIRST 2 OF EACH
%     'HERO_gka1'


session = 'all' ;
% *** Dates ***
%     '041416' ...
%     '041516' ...
        
% 1 -> use canonical HRF, 0 -> extract HRF using FIR
bCanonicalHRF = 0;

% Boolean: 1 -> go into debug mode--only fit light flux A
bDEBUG = 1;

bFreeFloatParams = 0;
%% LOAD TIME SERIES AND GET STIMULUS (& ATTENTION) START TIMES

% load time series
[avgTS, avgTSprc, tsFileNames, stimTypeArr, runOrder] ...
= loadTimeSeriesData(subj_name,session);

% get all stimulus values and start times, as well as the attention task
% start times
[startTimesSorted, stimValuesSorted, attnStartTimes] = orderStartTimes(subj_name,session);

% Time Series sampling points
timebase = [1:336]-1;

%% HRF PARAMETERS

% how long we expect the HRF to be
lengthHRF = 15;

%% DERIVE HRF FROM DATA, CREATE STIMULUS MODELS

% derive HRF from data
[BOLDHRF, cleanedData, SEHRF]= deriveHRFwrapper(avgTSprc,attnStartTimes,lengthHRF,'Fourier');

% in case we use the FIR extracted HRF; if we are not, 'hrf' never gets
% used
hrfPointsToSample = [1 1000:1000:length(BOLDHRF)];

hrf = BOLDHRF(hrfPointsToSample);

if bCanonicalHRF == 1     
   % Double Gamma HRF--get rid of the FIR-extracted HRF from earlier
   clear BOLDHRF
   clear hrf
   BOLDHRF = createCanonicalHRF(0:lengthHRF,6,12,10);
else 
   % initialize vector for HRF
   BOLDHRF = zeros([1 size(avgTSprc,2)]);
   % align HRF with 0 mark
   hrf = hrf-hrf(1);
   
   figure;
   errorbar(0:lengthHRF,hrf,SEHRF(hrfPointsToSample),'LineWidth',2)
   xlabel('Time/s'); ylabel('Signal'); set(gca,'FontSize',15);
   title('HRF');
   
   % make it the right size
   BOLDHRF(1:length(hrf)) = hrf;     
end

%% STIMULUS VECTOR CREATION

% resolution to sample stimulus step function
stepFunctionRes = 50;
% length of cosine ramp (seconds)
cosRamp = 3;
% stimulus duration
stimDuration = 12;

if bFreeFloatParams == 1
   startTimesSorted = 0:stimDuration:max(timebase); 
   startTimesSorted = repmat(startTimesSorted,[size(stimValuesSorted,1) 1]);
end

%% GET BETA AND MODEL FIT

% create stimulus vector
[stimMatrix,stimValues,startTimesSorted_A,startTimesSorted_B, ...
stimValuesSorted_A,stimValuesSorted_B,actualStimulusValues] ...
= createStimMatrix(startTimesSorted,stimValuesSorted,tsFileNames, ...
timebase,stimDuration,stepFunctionRes,cosRamp,bFreeFloatParams);

if bFreeFloatParams == 1
    stimValuesSorted_A = [0 stimValuesSorted_A];    
end

%% PRELIMINARIES FOR FITTING

% store parameters--initialize matrices
storeUnique = zeros([0 0 0]);
storeAll = zeros([0 0 0]);

reconstructedTSmat = [];
MSEstore = [];

% if in debug mode, only fit light flux A
if bDEBUG == 1
   runsToFit = find(stimTypeArr == 1 & runOrder == 'A');
else
   runsToFit = 1:size(stimMatrix,1);  
end

for i = 1:length(runsToFit)
   % wrap various things into a stimulus struct
   stimulus = struct;
   stimulus.stimMatrix = squeeze(stimMatrix(runsToFit(i),:,:));
   stimulus.stimValues = squeeze(stimValues(runsToFit(i),:));
   stimulus.startTimesSorted_A = startTimesSorted_A;
   stimulus.startTimesSorted_B = startTimesSorted_B;
   stimulus.stimValuesSorted_A = stimValuesSorted_A;
   stimulus.stimValuesSorted_B = stimValuesSorted_B;
   stimulus.uniqueTemporalFreq = unique(stimValuesSorted_A);
   % pull out one BOLD response
   boldResponse = squeeze(cleanedData(runsToFit(i),:));
   % number of stimuli
   defaultParamsInfo.nStimuli = size(stimulus.stimMatrix,1);
   % initialize object
   tmri = tmriTFBlockDesignColorModel;
   % set initial parameters
   params0 = tmri.defaultParams('DefaultParamsInfo',defaultParamsInfo);
   fprintf('Default model parameters:\n');
   tmri.print(params0);
   % HRF struct
   theHRF.HRF = BOLDHRF;
   theHRF.timebase = timebase;
   % get parameter locking matrix
   paramLockMatrix = tmri.lockMatrix(params0,stimulus);
   % fit the response
   [paramsFit,fitResponse,fval] = tmri.fitResponse(timebase,stimulus,boldResponse, ...
                          'HRF',theHRF,'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
   fprintf('Model parameter from fits:\n');
   tmri.print(paramsFit);
   % store all parameters for each run, regardless of daisy-chaining
   storeAll(size(storeAll,1)+1,:,:) = paramsFit.paramMainMatrix;     
   % store MSE measure
   MSEstore(length(MSEstore)+1) = fval;
   [~,meanParamValues] = tmri.plotParams(paramsFit,stimulus);
   % do this for each run
   storeUnique(size(storeUnique,1)+1,:,:) = meanParamValues;
   reconstructedTSmat(size(reconstructedTSmat,1)+1,:) = fitResponse;
   display(['run number: ' num2str(i)]);
end
