% This Script creates Temporal Transfer Functions 
% using a linear model and FIR with cosine windows.

%% Variable name legend 

% ts        = time series
% LH & RH   = left & right hemisphere
% stim      = stimulus
% AVG       = average

%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_gka1' ; 
% *** Subject Pool ***
%     'HERO_asb1' 
%     'HERO_gka1'


session = 'all' ;
% *** Dates ***
%     '041416' ...
%     '041516' ...
        
bCanonicalHRF = 0;

%% Defining Paths, Order, etc

% load time series
[avgTS, avgTSprc, tsFileNames, stimTypeArr, runOrder] ...
= loadTimeSeriesData(subj_name,session);

% get all stimulus values and start times, as well as the attention task
% start times
[startTimesSorted, stimValuesSorted, attnStartTimes] = orderStartTimes(subj_name,session);

% how long we expect the HRF to be
lengthHRF = 26;
% acquisition time
T_R = 1;
% resolution to sample stimulus step function
stepFunctionRes = 50;
% length of cosine ramp (seconds)
cosRamp = 3;
% Time Series sampling points
TS_timeSamples = 1:336;
% stimulus duration
stimDuration = 12;

% derive HRF from data
[BOLDHRF, cleanedData]= deriveHRF(avgTSprc,attnStartTimes,lengthHRF,T_R);

% get unique stimulus values
actualStimulusValues = unique(stimValuesSorted(stimValuesSorted~=-1 & stimValuesSorted~=0));

% Stire Stimulus Order A & B
stimValuesSorted_A = [] ;
stimValuesSorted_B = [] ;

% for each run
for i = 1:size(startTimesSorted,1)
    % Stores A & B sequences-- For labeling Time Series by Stimulus Period
   if strfind(char(tsFileNames(i)),'_A_') & isempty(stimValuesSorted_A)
      startTimesSorted_A = startTimesSorted(i,startTimesSorted(i,:)~=-1);  
      stimValuesSorted_A = stimValuesSorted(i,stimValuesSorted(i,:)~=-1); 
   elseif strfind(char(tsFileNames(i)),'_B_') & isempty(stimValuesSorted_B)
      startTimesSorted_B = startTimesSorted(i,startTimesSorted(i,:)~=-1); 
      stimValuesSorted_B = stimValuesSorted(i,stimValuesSorted(i,:)~=-1); 
   else
      stimOrderMarker = [] ; 
   end
    % and each stimulus
   for j = 1:length(actualStimulusValues)
       % grab the starting times
       startTimesForGivenStimValue = startTimesSorted(i,stimValuesSorted(i,:)==actualStimulusValues(j));
       % create stimulus model for each one, and sum those models together
       % to get the stimulus model for each stimulus type
       for k = 1:length(startTimesForGivenStimValue)
           halfCosine(k,:) = createStimVector(TS_timeSamples,startTimesForGivenStimValue(k), ...
                        stimDuration,stepFunctionRes,cosRamp);
       end
       vectorForOneStim(i,j,:) = sum(halfCosine);
   end
end

% Total Duration is simply Largest time value
modelDuration=floor(max(TS_timeSamples)) ; 
modelResolution=20 ; 

% Time Samples to Interpolate
t = linspace(1,modelDuration,modelDuration.*modelResolution) ;

% in case we use the FIR extracted HRF 
if strcmp(subj_name,'HERO_asb1')
  hrf = BOLDHRF(1:10);
elseif strcmp(subj_name,'HERO_gka1')
  hrf = BOLDHRF(1:14);
else
  error('BOLDmodelFitScript: invalid subject');
end
       
if bCanonicalHRF == 1     
   % Double Gamma HRF
   clear BOLDHRF
   BOLDHRF = createCanonicalHRF(t,6,12,10);
else 
   BOLDHRF_unInterp = zeros([1 length(avgTS(i,:))]);
   % ALIGN HRF WITH 0 MARK
   hrf = hrf-hrf(1);
   BOLDHRF_unInterp(1:length(hrf)) = hrf;
   % UPSAMPLE HRF
   BOLDHRF = interp1(TS_timeSamples,BOLDHRF_unInterp,t);
   BOLDHRF(isnan(BOLDHRF)) = 0;       
end

% LOOP OVER RUNS
for i = 1:size(vectorForOneStim,1)
    % LOOP OVER STIMULI WITHIN EACH RUN
   for j = 1:size(vectorForOneStim,2)
       regressor = createRegressor(squeeze(vectorForOneStim(i,j,:))',TS_timeSamples,BOLDHRF,t);
       designMatrixPreOnes(:,j) = regressor-mean(regressor);
   end
   designMatrix = [ones([size(designMatrixPreOnes,1) 1]) designMatrixPreOnes] ;
    % Obtain Beta Weights & Plot
   betaWeights = designMatrix\cleanedData(i,:)' ; 

   % Beta Weights Sans Weight for the first Regressor
   betaMatrix(i,:) = betaWeights(2:length(betaWeights)) ;

   % Reconstruct Time Series According to Model 
   reconstructedTS = sum(repmat(betaWeights',[size(designMatrix,1) 1]).*designMatrix,2) ;

   % Store all reconstructed Time Series
   reconstructedTSmat(i,:) = reconstructedTS ;
end

% Self-Explanatory Variable Names
numberOfRuns = 12 ;
numRunsPerStimOrder = 6 ;   % Stim order A -or- B

%% Calculations (Means and Standard Errors)

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
