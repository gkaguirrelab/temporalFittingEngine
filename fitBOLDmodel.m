% This Script creates Temporal Transfer Functions 
% using a linear model and FIR with cosine windows.

%% Variable name legend 

% ts        = time series
% LH & RH   = left & right hemisphere
% stim      = stimulus
% AVG       = average

%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_asb1' ; 
% *** Subject Pool ***
%     'HERO_asb1' 
%     'HERO_gka1'


session = 'all' ;
% *** Dates ***
%     '041416' ...
%     '041516' ...
        
% 1 -> use canonical HRF, 0 -> extract HRF using FIR
bCanonicalHRF = 0;

%% LOAD TIME SERIES AND GET STIMULUS (& ATTENTION) START TIMES

% load time series
[avgTS, avgTSprc, tsFileNames, stimTypeArr, runOrder] ...
= loadTimeSeriesData(subj_name,session);

% get all stimulus values and start times, as well as the attention task
% start times
[startTimesSorted, stimValuesSorted, attnStartTimes] = orderStartTimes(subj_name,session);

% Time Series sampling points
TS_timeSamples = 1:336;

%% HRF PARAMETERS

% how long we expect the HRF to be
lengthHRF = 26;

% acquisition time (only useful for HRF so far)
T_R = 1;

% These parameters pertain to the HRF. Total Duration is simply Largest 
% time value
modelDuration=floor(max(TS_timeSamples)) ; 
modelResolution=20 ; 

% Time Samples to Interpolate
t = linspace(1,modelDuration,modelDuration.*modelResolution) ;

%% DERIVE HRF FROM DATA, CREATE STIMULUS MODELS

% derive HRF from data
[BOLDHRF, cleanedData]= deriveHRF(avgTSprc,attnStartTimes,lengthHRF,T_R);

% in case we use the FIR extracted HRF; if we are not, 'hrf' never gets
% used
if strcmp(subj_name,'HERO_asb1')
  hrf = BOLDHRF(1:10);
elseif strcmp(subj_name,'HERO_gka1')
  hrf = BOLDHRF(1:14);
else
  error('BOLDmodelFitScript: invalid subject');
end
       
if bCanonicalHRF == 1     
   % Double Gamma HRF--get rid of the FIR-extracted HRF from earlier
   clear BOLDHRF
   clear hrf
   BOLDHRF = createCanonicalHRF(t,6,12,10);
else 
   % initialize vector for HRF
   BOLDHRF_unInterp = zeros([1 size(avgTSprc,2)]);
   % align HRF with 0 mark
   hrf = hrf-hrf(1);
   % make it the right size
   BOLDHRF_unInterp(1:length(hrf)) = hrf;
   % upsample the HRF
   BOLDHRF = interp1(TS_timeSamples,BOLDHRF_unInterp,t);
   BOLDHRF(isnan(BOLDHRF)) = 0;       
end

%% STIMULUS VECTOR CREATION

% resolution to sample stimulus step function
stepFunctionRes = 50;
% length of cosine ramp (seconds)
cosRamp = 3;
% stimulus duration
stimDuration = 12;

%% GET BETA AND MODEL FIT

[betaMatrix,reconstructedTSmat,startTimesSorted_A, ...
          stimValuesSorted_A, startTimesSorted_B, stimValuesSorted_B, actualStimulusValues] ...
                                      = ...
          forwardModel(startTimesSorted,stimValuesSorted,tsFileNames, ...
          TS_timeSamples,stimDuration,stepFunctionRes,cosRamp, ...
          t,BOLDHRF,cleanedData);

%% Calculations (Means and Standard Errors)

% Self-Explanatory Variable Names
numberOfRuns = 12 ;
numRunsPerStimOrder = 6 ;   % Stim order A -or- B

% Convert Mean-Subtracted Beta values to Percentages
LightFluxBeta =  mean(betaMatrix(stimTypeArr == 1,:));
L_minus_M_Beta = mean(betaMatrix(stimTypeArr == 2,:));
S_Beta =         mean(betaMatrix(stimTypeArr == 3,:));

% Compute Standard Error
LightFluxBetaSE =  ((std(betaMatrix(stimTypeArr == 1,:)))./sqrt(numberOfRuns));
L_minus_M_BetaSE = ((std(betaMatrix(stimTypeArr == 2,:)))./sqrt(numberOfRuns));
S_BetaSE =         ((std(betaMatrix(stimTypeArr == 3,:)))./sqrt(numberOfRuns));

% Average Time Series for Each Combination of Stimulus Type & Run order
LightFluxAvgTS_A =  mean(cleanedData(stimTypeArr == 1 & runOrder == 'A',:)) ;
L_minus_M_AvgTS_A = mean(cleanedData(stimTypeArr == 2 & runOrder == 'A',:)) ;
S_AvgTS_A =         mean(cleanedData(stimTypeArr == 3 & runOrder == 'A',:)) ;

LightFluxAvgTS_B =  mean(cleanedData(stimTypeArr == 1 & runOrder == 'B',:)) ;
L_minus_M_AvgTS_B = mean(cleanedData(stimTypeArr == 2 & runOrder == 'B',:)) ;
S_AvgTS_B =         mean(cleanedData(stimTypeArr == 3 & runOrder == 'B',:)) ;

% Standard Error of Time Series
LightFluxStdTS_A =  (std(cleanedData(stimTypeArr == 1 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;
L_minus_M_StdTS_A = (std(cleanedData(stimTypeArr == 2 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;
S_StdTS_A =         (std(cleanedData(stimTypeArr == 3 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;

LightFluxStdTS_B =  (std(cleanedData(stimTypeArr == 1 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;
L_minus_M_StdTS_B = (std(cleanedData(stimTypeArr == 2 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;
S_StdTS_B =         (std(cleanedData(stimTypeArr == 3 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;

% Do the Same for 'Reconstructed' Time Series
LightFluxAvgTS_Model_A =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'A',:)) ;
L_minus_M_AvgTS_Model_A = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'A',:)) ;
S_AvgTS_Model_A =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'A',:)) ;

LightFluxAvgTS_Model_B =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'B',:)) ;
L_minus_M_AvgTS_Model_B = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'B',:)) ;
S_AvgTS_Model_B =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'B',:)) ;

yLimits = [min([LightFluxBeta L_minus_M_Beta S_Beta]) max([LightFluxBeta L_minus_M_Beta S_Beta])] ;

%% TTF & HRF Plots

% Light Flux
[wftd1, fp1] = fitWatsonToTTF_errorGuided(actualStimulusValues',LightFluxBeta,LightFluxBetaSE,1); hold on
errorbar(actualStimulusValues',LightFluxBeta,LightFluxBetaSE,'ko'); set(gca,'FontSize',15);
set(gca,'Xtick',actualStimulusValues'); title('Light flux');

[wftd2, fp2] = fitWatsonToTTF_errorGuided(actualStimulusValues',L_minus_M_Beta,L_minus_M_BetaSE,1); hold on
errorbar(actualStimulusValues',L_minus_M_Beta,L_minus_M_BetaSE,'ko');
set(gca,'FontSize',15); set(gca,'Xtick',actualStimulusValues'); title('L - M');
% S
[wftd3, fp3] = fitWatsonToTTF_errorGuided(actualStimulusValues',S_Beta,S_BetaSE,1); hold on
errorbar(actualStimulusValues',S_Beta,S_BetaSE,'ko'); set(gca,'FontSize',15);
set(gca,'Xtick',actualStimulusValues'); title('S');

%% Time Series plots 
% Use Function for plotting Data:
% -- plotLinModelFits -- 

% Create Cells for Labeling Plots
stimValuesMatSorted_A_cell = {} ;
for i = 1:length(stimValuesSorted_A)
   stimValuesMatSorted_A_cell{i} = num2str(stimValuesSorted_A(i)) ; 
end

stimValuesMatSorted_B_cell = {} ;
for i = 1:length(stimValuesSorted_B)
   stimValuesMatSorted_B_cell{i} = num2str(stimValuesSorted_B(i)) ; 
end

% Set Figure Dimensions
figure;
set(gcf,'Position',[156 372 1522 641])

% Light Flux -A
subplot(3,2,1)
plotLinModelFits(T_R.*(1:length(LightFluxAvgTS_A)),LightFluxAvgTS_A,LightFluxAvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,LightFluxStdTS_A);
title('Light flux A'); xlabel('Time / s'); ylabel('% signal change');

% L minus M -A
subplot(3,2,3)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_A)),L_minus_M_AvgTS_A,L_minus_M_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,L_minus_M_StdTS_A);
title('L - M A'); xlabel('Time / s'); ylabel('% signal change');

% S -A
subplot(3,2,5)
plotLinModelFits(T_R.*(1:length(S_AvgTS_A)),S_AvgTS_A,S_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,S_StdTS_A);
title('S A'); xlabel('Time / s'); ylabel('% signal change');

% Light Flux -B
subplot(3,2,2)
plotLinModelFits(T_R.*(1:length(LightFluxAvgTS_B)),LightFluxAvgTS_B,LightFluxAvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,LightFluxStdTS_B);
title('Light flux B');

% L minus M -B
subplot(3,2,4)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_B)),L_minus_M_AvgTS_B,L_minus_M_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,L_minus_M_StdTS_B);
title('L - M B');

% S -B
subplot(3,2,6)
plotLinModelFits(T_R.*(1:length(S_AvgTS_B)),S_AvgTS_B,S_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,S_StdTS_B);
title('S B');
