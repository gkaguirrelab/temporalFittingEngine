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
lengthHRF = 16;

% acquisition time (only useful for HRF so far)
T_R = 1;

% These parameters pertain to the HRF. Total Duration is simply Largest 
% time value
modelDuration=floor(max(TS_timeSamples)) ; 
modelSampleFreq=20 ; 

% Time Samples to Interpolate
modelUpsampled_t = linspace(1,modelDuration,modelDuration.*modelSampleFreq) ;

%% DERIVE HRF FROM DATA, CREATE STIMULUS MODELS

% derive HRF from data
[BOLDHRF, cleanedData]= fitHRF_FIR(avgTSprc,attnStartTimes,lengthHRF,T_R);

% in case we use the FIR extracted HRF; if we are not, 'hrf' never gets
% used
if strcmp(subj_name,'HERO_asb1')
  hrf = BOLDHRF(1:lengthHRF);
elseif strcmp(subj_name,'HERO_gka1')
  hrf = BOLDHRF(1:lengthHRF);
else
  error('BOLDmodelFitScript: invalid subject');
end
       
if bCanonicalHRF == 1     
   % Double Gamma HRF--get rid of the FIR-extracted HRF from earlier
   clear BOLDHRF
   clear hrf
   BOLDHRF = createCanonicalHRF(modelUpsampled_t,6,12,10);
else 
   % initialize vector for HRF
   BOLDHRF_unInterp = zeros([1 size(avgTSprc,2)]);
   % align HRF with 0 mark
   hrf = hrf-hrf(1);
   % make it the right size
   BOLDHRF_unInterp(1:length(hrf)) = hrf;
   % upsample the HRF
   BOLDHRF = interp1(TS_timeSamples,BOLDHRF_unInterp,modelUpsampled_t);
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

% create stimulus vector
[stimMatrix,paramLockMatrix,startTimesSorted_A,startTimesSorted_B, ...
stimValuesSorted_A,stimValuesSorted_B,actualStimulusValues] ...
= createStimMatrix(startTimesSorted,stimValuesSorted,tsFileNames, ...
TS_timeSamples,stimDuration,stepFunctionRes,cosRamp);

%%
% store the HRF, its time samples, and the amplitudes in the struct
paramStruct.HRF = BOLDHRF;
paramStruct.HRFtimeSamples = modelUpsampled_t;

paramStruct.Amplitude = 0.5.*ones([size(stimMatrix,2) 1]);
paramStruct.tau2 = 0.001.*ones([size(stimMatrix,2) 1]);
paramStruct.ARAmplitude = (-0.125).*ones([size(stimMatrix,2) 1]);

%%
% store amplitudes
ampStore = [];
tau2store = [];
ARampStore = [];

reconstructedTSmat = [];
MSEstore = [];

% Boolean: 1 -> go into debug mode--only fit light flux A
bDEBUG = 0;

if bDEBUG == 1
   runsToFit = find(stimTypeArr == 1 & runOrder == 'A');
else
   runsToFit = 1:size(stimMatrix,1);  
end

for i = 1:length(runsToFit)
    % call fitting routine
    [paramStructFit,fval]= fitNeuralParams(squeeze(stimMatrix(runsToFit(i),:,:)),TS_timeSamples,squeeze(paramLockMatrix(runsToFit(i),:,:)),cleanedData(runsToFit(i),:),paramStruct);
    amp = paramStructFit.Amplitude;
    ARamp = paramStructFit.ARAmplitude;
    tau2forStim = paramStructFit.tau2;
    MSEstore(length(MSEstore)+1) = fval;
     % Determine which stimulus values went with which parameter
   if strfind(char(tsFileNames(runsToFit(i))),'_A_')
      valueLookup = stimValuesSorted_A(stimValuesSorted_A>0);
   elseif strfind(char(tsFileNames(runsToFit(i))),'_B_')
      valueLookup = stimValuesSorted_B(stimValuesSorted_B>0);
   else
      valueLookup = [] ; 
   end
   % get only unique stim values, and their corresponding locked params
   [stimValueToPlot,ia] = unique(valueLookup);
   
    % store fit amplitudes 
    ampStore(size(ampStore,1)+1,:) = amp(ia);
    tau2store(size(tau2store,1)+1,:) = tau2forStim(ia);
    ARampStore(size(ARampStore,1)+1,:) = ARamp(ia);
    % store reconstructed time series
     [~,reconstructedTS] = forwardModel(squeeze(stimMatrix(runsToFit(i),:,:)),TS_timeSamples,cleanedData(runsToFit(i),:),paramStructFit);
     reconstructedTSmat(size(reconstructedTSmat,1)+1,:) = reconstructedTS;
end
%%
if bDEBUG == 1
   Beta = median(ampStore); 
   BetaSE = std(ampStore)./sqrt(size(ampStore,1));
   tau2 = median(tau2store); 
   tau2SE = std(tau2store)./sqrt(size(tau2store,1));
   AR = median(ARampStore); 
   ARSE = std(ARampStore)./sqrt(size(ARampStore,1));
   AvgTS = mean(cleanedData(runsToFit,:));
   StdTS = std(cleanedData(runsToFit,:))./sqrt(size(ampStore,1));
   MSE = mean(MSEstore);
   AvgTS_model = mean(reconstructedTSmat);
   
   stimValuesMatSorted_A_cell = {} ;
    for j = 1:length(stimValuesSorted_A)
       stimValuesMatSorted_A_cell{j} = num2str(stimValuesSorted_A(j)) ; 
    end
   
   figure;
   errorbar(actualStimulusValues',Beta,BetaSE,'-ko'); set(gca,'FontSize',15);
    set(gca,'Xtick',actualStimulusValues'); title('Light flux'); set(gca,'Xscale','log');
    xlabel('Temporal frequency (Hz)'); ylabel('median main-response amplitude');
    
    figure;
    errorbar(actualStimulusValues',tau2,tau2SE,'-ko'); set(gca,'FontSize',15);
    set(gca,'Xtick',actualStimulusValues'); title('Light flux'); set(gca,'Xscale','log');
    xlabel('Temporal frequency (Hz)'); ylabel('median \tau_2');

    figure;
    errorbar(actualStimulusValues',AR,ARSE,'-ko'); set(gca,'FontSize',15);
    set(gca,'Xtick',actualStimulusValues'); title('Light flux'); set(gca,'Xscale','log');
    xlabel('Temporal frequency (Hz)'); ylabel('median after-response amplitude');    

    figure;
    plotLinModelFits(T_R.*(1:length(AvgTS)),AvgTS,AvgTS_model, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,StdTS,MSE);
    title('Light flux A'); xlabel('Time / s'); ylabel('% signal change');
end

%%
% Self-Explanatory Variable Names
numberOfRuns = 12 ;
numRunsPerStimOrder = 6 ;   % Stim order A -or- B

%% Calculations (Means and Standard Errors)

% Convert Mean-Subtracted Beta values to Percentages
LightFluxBeta =  mean(ampStore(stimTypeArr == 1,:));
L_minus_M_Beta = mean(ampStore(stimTypeArr == 2,:));
S_Beta =         mean(ampStore(stimTypeArr == 3,:));

% Compute Standard Error
LightFluxBetaSE =  ((std(ampStore(stimTypeArr == 1,:)))./sqrt(numberOfRuns));
L_minus_M_BetaSE = ((std(ampStore(stimTypeArr == 2,:)))./sqrt(numberOfRuns));
S_BetaSE =         ((std(ampStore(stimTypeArr == 3,:)))./sqrt(numberOfRuns));

% Convert Mean-Subtracted Beta values to Percentages
LightFluxtau2 =  mean(tau2store(stimTypeArr == 1,:));
L_minus_M_tau2 = mean(tau2store(stimTypeArr == 2,:));
S_tau2 =         mean(tau2store(stimTypeArr == 3,:));

% Compute Standard Error
LightFluxtau2SE =  ((std(tau2store(stimTypeArr == 1,:)))./sqrt(numberOfRuns));
L_minus_M_tau2SE = ((std(tau2store(stimTypeArr == 2,:)))./sqrt(numberOfRuns));
S_tau2SE =         ((std(tau2store(stimTypeArr == 3,:)))./sqrt(numberOfRuns));

%%
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
%% MEAN SQUARED ERROR VALUES
LightFluxMSE_A =  mean(MSEstore(stimTypeArr == 1 & runOrder == 'A')) ;
L_minus_M_MSE_A = mean(MSEstore(stimTypeArr == 2 & runOrder == 'A')) ;
S_MSE_A =         mean(MSEstore(stimTypeArr == 3 & runOrder == 'A')) ;

LightFluxMSE_B =  mean(MSEstore(stimTypeArr == 1 & runOrder == 'B')) ;
L_minus_M_MSE_B = mean(MSEstore(stimTypeArr == 2 & runOrder == 'B')) ;
S_MSE_B =         mean(MSEstore(stimTypeArr == 3 & runOrder == 'B')) ;

%%

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

%% TAU VALUES
figure;
errorbar(actualStimulusValues,LightFluxtau2,LightFluxtau2SE,'ko'); set(gca,'FontSize',15); hold on
set(gca,'Xtick',actualStimulusValues'); title('Light flux');
errorbar(actualStimulusValues,L_minus_M_tau2,L_minus_M_tau2SE,'ro'); set(gca,'FontSize',15);
errorbar(actualStimulusValues,S_tau2,S_tau2SE,'bo'); set(gca,'FontSize',15);
set(gca,'Xscale','log');
xlabel('Temporal frequency'); ylabel('\tau_2');

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
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,LightFluxStdTS_A,LightFluxMSE_A);
title('Light flux A'); xlabel('Time / s'); ylabel('% signal change');

% L minus M -A
subplot(3,2,3)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_A)),L_minus_M_AvgTS_A,L_minus_M_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,L_minus_M_StdTS_A,L_minus_M_MSE_A);
title('L - M A'); xlabel('Time / s'); ylabel('% signal change');

% S -A
subplot(3,2,5)
plotLinModelFits(T_R.*(1:length(S_AvgTS_A)),S_AvgTS_A,S_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,S_StdTS_A,S_MSE_A);
title('S A'); xlabel('Time / s'); ylabel('% signal change');

% Light Flux -B
subplot(3,2,2)
plotLinModelFits(T_R.*(1:length(LightFluxAvgTS_B)),LightFluxAvgTS_B,LightFluxAvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,LightFluxStdTS_B,LightFluxMSE_B);
title('Light flux B');

% L minus M -B
subplot(3,2,4)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_B)),L_minus_M_AvgTS_B,L_minus_M_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,L_minus_M_StdTS_B,L_minus_M_MSE_B);
title('L - M B');

% S -B
subplot(3,2,6)
plotLinModelFits(T_R.*(1:length(S_AvgTS_B)),S_AvgTS_B,S_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,S_StdTS_B,S_MSE_B);
title('S B');
