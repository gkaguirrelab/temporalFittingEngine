function [LightFluxAvgTS_A, L_minus_M_AvgTS_A, S_AvgTS_A, ...
          LightFluxAvgTS_B, L_minus_M_AvgTS_B, S_AvgTS_B, ...
          LightFluxStdTS_A, L_minus_M_StdTS_A, S_StdTS_A, ...
          LightFluxStdTS_B, L_minus_M_StdTS_B, S_StdTS_B, ...
          LightFluxMSE_A, L_minus_M_MSE_A, S_MSE_A, ...
          LightFluxMSE_B, L_minus_M_MSE_B, S_MSE_B, ...
          LightFluxAvgTS_Model_A, L_minus_M_AvgTS_Model_A, S_AvgTS_Model_A, ...
          LightFluxAvgTS_Model_B, L_minus_M_AvgTS_Model_B, S_AvgTS_Model_B] ...
         = ...
summaryStatsOfFits(cleanedData,MSEstore,reconstructedTSmat,stimTypeArr,runOrder)

% creates summary statistics for fit plots for BDCM. ***SPECIFIC TO
% BDCM***, i.e. not generalizable. mainly for making fitBOLDmodel more
% presentable

% store variables
numberOfRuns = sum(stimTypeArr==1) ;
numRunsPerStimOrder = sum(stimTypeArr==1 & runOrder=='A') ;   % Stim order A -or- B    

% TIME SERIES MEAN AND STD ERROR

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

%% MEANS FOR MODEL FITS

% Do the Same for 'Reconstructed' Time Series
LightFluxAvgTS_Model_A =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'A',:)) ;
L_minus_M_AvgTS_Model_A = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'A',:)) ;
S_AvgTS_Model_A =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'A',:)) ;

LightFluxAvgTS_Model_B =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'B',:)) ;
L_minus_M_AvgTS_Model_B = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'B',:)) ;
S_AvgTS_Model_B =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'B',:)) ;
    
end