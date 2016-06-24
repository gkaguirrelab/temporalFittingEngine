%% temporalTransferFunc

% This Script creates Temporal Transer Functions 
% using a linear model and FIR with cosine windows.

%% Variable name legend 

% ts        = time series
% LH & RH   = left & right hemisphere
% stim      = stimulus
% AVG       = average
% Arr       = array (superfluous--one of Ben's coding habit)
% Mat       = matrix
% cur       = current
% position  = index

%% Identify the user
 if isunix
    [~, user_name] = system('whoami') ; % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%') ; % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
 end


%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_asb1' ; 
% *** Subject Pool ***
%     'HERO_asb1' 
%     'HERO_gka1'


session = 'all' ;
% *** Dates ***
%     '041416' ...
%     '041516' ...

% Path to local Dropbox
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'] ;

% Define path to folder for one sunject on one date
dirPathStim = [localDropboxDir 'MELA_analysis/HCLV_Photo_7T/mriTemporalFitting_data/' ...
           subj_name '/' session '/' 'Stimuli/'] ;

%% HRF Parameters --(Grabbed from Winawer Model Code)

bCanonicalHRF = 1 ;
% 1 = Use Canonical HRF
% 0 = Use FIR-extracted HRF

if bCanonicalHRF == 1
    param = struct ;
    % parameters of the double-gamma hemodynamic filter (HRF)
    param.gamma1 = 6 ;      % positive gamma parameter (roughly, time-to-peak in seconds)
    param.gamma2 = 12 ;     % negative gamma parameter (roughly, time-to-peak in seconds)
    param.gammaScale = 10 ; % scaling factor between the positive and negative gamma componenets
end


% Duration of Stimulus (Constant)
stimTime = 12 ;

% Duration of Attention Task (Constant)
attnTime = 0.25 ;

% Attention Task has no Frequency, but it gets a 'code'
attnCode = 1 ;

% Boolean to determine if we want to model the Attention Task with FIR
attnFIR = 1 ;  % 1 = YES :: 0 = NO

% Length ofAttention HRF
lengthAttnHRF = 26 ;

% Acquisition Time
T_R = 1 ;

% Time Series sampling points
TS_timeSamples = 1:336 ;
        
%% Defining Paths, Order, etc

% load time series
[avgTS, timeSeriesStore, tsFileNames, stimTypeArr, runOrder] = loadTimeSeriesData(subj_name,session);

%% Load & Plot Stimulus Step Functions

% Load all contents of Stimulus Directory
stimDirContents = dir(dirPathStim) ;

% Number of Stimulus Folders
numberOfFolders = length(stimDirContents) ;

% Initialize Cell containing all Stimulus folder names
folderNameCell = {} ;

% Loop over numebr of Stimulus folders & create cell with their names
for i = 1:numberOfFolders
   miniFolderName = stimDirContents(i).name ;
   if length(miniFolderName)>4 & strcmp(miniFolderName(1:4),'HERO') ;
       folderNameCell{length(folderNameCell)+1} = miniFolderName ;
   end
end

% Initialize Matrix for storing Beta values
betaMatrix = [] ;

% Stire Stimulus Order A & B
stimValuesSorted_A = [] ;
stimValuesSorted_B = [] ;

% Loop Over Stimulus folder names
for i = 1:length(folderNameCell)
   % Look in each run's folder
   currentDirPath = [dirPathStim char(folderNameCell(i))] ; 
   % Get all their contents
   runFiles = dir(currentDirPath) ;

   % Initialize Matrices for storing start Times & Stimulus values
   startTimes = [] ;
   stimValues = [] ;
   
   attnFunctionRes = 10 ; % Attention Function Resolution
   
   % Loop over files in each folder
   for j = 1:length(runFiles)
       lengthOfCurFile = length(runFiles(j).name);
       curFile = runFiles(j).name;

       % We are interested in Hz_all files
       if length(curFile)>10 & strcmp(curFile(length(curFile)-9:length(curFile)),'Hz_all.txt')
          % Extract Temporal Frequency of Stimulus from file name
          stimFile = load([currentDirPath '/' curFile]) ; 
          freqValueTxt = curFile(length(curFile)-11:length(curFile)-10) ;

          % Number (0,2,4,8,16,32,64)-- can be 2-digit or 1-digit
          if str2num(freqValueTxt)
              freqValueNum = str2num(freqValueTxt) ;
          else
              freqValueNum = str2num(freqValueTxt(2)) ;            
          end

          % Grab all values in first column (Starting times)
          curTimeValue = stimFile(:,1) ;
          % Collect all start times and corresponding Stimulus values
          startTimes(length(startTimes)+1:length(startTimes)+length(curTimeValue)) = curTimeValue ;
          stimValues(length(stimValues)+1:length(stimValues)+length(curTimeValue)) = freqValueNum ;         
          
          % If the file contains Attention Task data
       elseif length(curFile)>20 & strcmp(curFile(length(curFile)-16:length(curFile)),'attentionTask.txt')
           attnFile = load([currentDirPath '/' curFile]) ;

           % Run FIR on Attention Start times
           [hrf,~] = attentionFIR(T_R.*(1:length(avgTS(i,:))),((avgTS(i,:)-mean(avgTS(i,:)))./mean(avgTS(i,:))).*100,attnFile(:,1),lengthAttnHRF,T_R) ;
           hrfStore(i,:) = hrf ;
           
           attnTimeValues = [] ;
           attnStimValues = [] ;
           
           % Anchor Attention Step function at 0
           attnStartTimes1 = attnFile(1,1) ;
           
           if abs(attnStartTimes1) > 0.01
              attnTimeValues(1) = 0 ;
              attnStimValues(1) = 0 ;
           end

           % Attention task is just a Dirac Delta Function (-ish)
           for k = 1:size(attnFile,1)
                curTimeValue = attnFile(k,1);
                attnTimeValues = [attnTimeValues [curTimeValue linspace(curTimeValue+1e-10,curTimeValue+attnTime-1e-3,attnFunctionRes) ...
                                 curTimeValue+attnTime-1e-7] ]; 
                attnStimValues = [attnStimValues [0 repmat(attnCode,[1 attnFunctionRes]) 0]] ;
           end
       end
       
   end

   % Sort the Big Vector of Start Times
   [startTimesSorted, stmsInd] = sort(startTimes) ;
   % Sort Corresponding Simulus values
   stimValuesSorted = stimValues(stmsInd) ;
   
   % Initialize Matrices for storing actual Stimulus plots
   stimPlotsTimeSamples = [] ;
   % Marks Which Stimulus is at each Time Point
   stimPlotsValues = [] ;
   % Stores the Cosine-Windowed Stimulus
   cosWindowedStim = [] ;
   
   % Make sure plot starts at 0 (So convolution won't return nans)
   if abs(startTimesSorted(1)) > 0.001 ;
      stimValuesSorted = [0 stimValuesSorted] ;
      startTimesSorted = [0 startTimesSorted] ;
   end

   % Stores A & B sequences-- For labeling Time Series by Stimulus Period
   if strfind(char(tsFileNames(i)),'_A_') & isempty(stimValuesSorted_A)
      startTimesSorted_A = startTimesSorted ;   
      stimValuesSorted_A = stimValuesSorted;
   elseif strfind(char(tsFileNames(i)),'_B_') & isempty(stimValuesSorted_B)
      startTimesSorted_B = startTimesSorted ;
      stimValuesSorted_B = stimValuesSorted ;
   else
      stimOrderMarker = [] ; 
   end

   % Sampling point between stimulus Onset & offset
   stepFunctionRes = 50 ;

   % Duration of Cosine Ramp (Seconds)
   cosRamp = 3 ;

   % Loop over Start Times
   for j = 1:length(startTimesSorted)
       
       % Grab each individual Start Time
       curTimeValue2 = startTimesSorted(j) ;
       curStimValue2 = stimValuesSorted(j) ;
       
       % Make 'BOX' Add tiny offset to make sure interpolation works well
       timeValues = [curTimeValue2 linspace(curTimeValue2+1e-7,curTimeValue2+stimTime-1e-5,stepFunctionRes) ...
                     curTimeValue2+stimTime-1e-7] ; 
       stimValues = [-1 repmat(curStimValue2,[1 length(timeValues)-2]) -1] ;

       % Shift Time Values so that Stimulus Onset becomes the 0
       timeForCos = timeValues - curTimeValue2 ;

       % Remove Stimulus Onset markers (for now)
       timeForCos = timeForCos(2:length(timeForCos)-1) ;
       
       % Define positions for left & right ramp
       leftRampInd = timeForCos <= cosRamp ;
       rightRampInd = timeForCos >= stimTime - cosRamp ;
       
       % Create half Cosine
       halfCosine = ones([1 length(timeForCos)]) ;
       halfCosine(leftRampInd) = fliplr((cos(linspace(0,pi,sum(leftRampInd)))+1)/2) ;
       halfCosine(rightRampInd) = (cos(linspace(0,pi,sum(rightRampInd)))+1)/2 ;

       % Put in points corresponding to Stimulus Start and End
       halfCosine = [-1 halfCosine -1] ;
       
       % Stick all the above in matrices defined before the loop
       stimPlotsTimeSamples = [stimPlotsTimeSamples timeValues] ;
       stimPlotsValues = [stimPlotsValues stimValues] ;
       cosWindowedStim = [cosWindowedStim halfCosine] ;
   end

   % All possible Stimuli (Sans Attention Task)
   stimHz = [2 4 8 16 32 64] ;

   % Initialize Matrix of Regressors
   regMatrix = [] ;
   
   % Make HRF (Taken from GKA's Winawer Model code)
   
       % Total Duration is simply Largest time value
       modelDuration=floor(max(stimPlotsTimeSamples)) ; 
       modelResolution=20 ; 

       % Time Samples to Interpolate
       t = linspace(1,modelDuration,modelDuration.*modelResolution) ;
       
   if bCanonicalHRF == 1     
       % Double Gamma HRF
       BOLDHRF = gampdf(t, param.gamma1, 1) - ...
       gampdf(t, param.gamma2, 1)/param.gammaScale ;
       % scale to unit sum to preserve amplitude of y following convolution
       BOLDHRF = BOLDHRF/sum(BOLDHRF) ;
   else

       load([subj_name '_HRF']);
       BOLDHRF_unInterp = zeros([1 length(avgTS(i,:))]);
       % ALIGN HRF WITH 0 MARK
       hrf = hrf-hrf(1);
       BOLDHRF_unInterp(1:length(hrf)) = hrf;
       % UPSAMPLE HRF
       BOLDHRF = interp1(TS_timeSamples,BOLDHRF_unInterp,t);
       BOLDHRF(isnan(BOLDHRF)) = 0;       
   end

   % Loop over Possible Stimulus values
   for j = 1:length(stimHz)
       
      % Get all positions with a given Stimulus 
      stimPositions = stimPlotsValues == stimHz(j) ;

      % Put in Cosine-Windowed Stimuli
      stimPreUps = zeros([1 length(stimPositions)]) ;
      stimPreUps(stimPositions) = cosWindowedStim(stimPositions) ;
      
      % Sample at Points t
      stimulusUpsampled = interp1(stimPlotsTimeSamples,stimPreUps,t,'linear','extrap') ;
      stimulusUpsampled = stimulusUpsampled(1:length(t)) ;
      
      % Convolve Stimulus with HRF to get Regressor
      regressorPreCut = conv(stimulusUpsampled,BOLDHRF) ;

      % Cut off extra Conv values --( Need to look more into this. Conv is
      % wierd in Matlab)
      regressor = regressorPreCut(1:length(stimulusUpsampled)) ;
      regressorDownSampled = interp1(t,regressor,TS_timeSamples,'linear','extrap') ;

      % Store Regressor
      regMatrix(:,j) = regressorDownSampled'-mean(regressorDownSampled) ; 
 
%        figure;
%        plot(stimPlotsTimeSamples,stimPositions); hold on
%        plot(t,stimulusUpsampled); 
%        xlabel('time/s');
%        plot(t,regressor); title(['Stimulus and BOLD signal for ' num2str(stimHz(j)) ' Hz' ' flicker']);
%        pause;
%        close;
   end
   
   % Add a Covariate of ones to the end
   regMatrix = [ones([size(regMatrix,1) 1]) regMatrix] ;

   % Get Step Function for Attention Task
   attnPositions = attnStimValues == attnCode ;
   attnPositions = double(attnPositions) ;

   % Sample it Evenly
   attnBox = interp1(attnTimeValues,attnStimValues,t) ;

   % Remove nans
   attnBox(isnan(attnBox)) = 0 ;
   
   % Convolve with HRF
   attnCovariate = conv(attnBox,BOLDHRF) ;
   
   % Sample to be the same length as other regressors
   attnCovariate = attnCovariate(1:length(t)) ;
   attnCovariateDownSampled = interp1(t,attnCovariate,TS_timeSamples,'linear','extrap') ;

   % Add to the Design Matrix
   regMatrix(:,size(regMatrix,2)+1) = attnCovariateDownSampled - mean(attnCovariateDownSampled) ;

   % Obtain Beta Weights & Plot
   betaWeights = regMatrix\avgTS(i,:)' ; 

   % Beta Weights Sans Weight for the first Regressor
   betaMatrix(i,:) = betaWeights(2:length(betaWeights)-1)./mean(avgTS(i,:)) ;

   % Reconstruct Time Series According to Model 
   % -- Convert to % signal change
   reconstructedTS = sum(repmat(betaWeights',[size(regMatrix,1) 1]).*regMatrix,2) ;
   reconstructedTS = ((reconstructedTS - mean(reconstructedTS))./mean(reconstructedTS)).*100 ;

   % Store all reconstructed Time Series
   reconstructedTSmat(i,:) = reconstructedTS ;

end

% Self-Explanatory Variable Names
numberOfRuns = 12 ;
numRunsPerStimOrder = 6 ;   % Stim order A -or- B

% Create Cells for Labeling Plots
stimValuesMatSorted_A_cell = {} ;
for i = 1:length(stimValuesSorted_A)
   stimValuesMatSorted_A_cell{i} = num2str(stimValuesSorted_A(i)) ; 
end

stimValuesMatSorted_B_cell = {} ;
for i = 1:length(stimValuesSorted_B)
   stimValuesMatSorted_B_cell{i} = num2str(stimValuesSorted_B(i)) ; 
end

%% Calculations (Means and Standard Errors)

% Convert Mean-Subtracted Beta values to Percentages
LightFluxBeta =  mean(betaMatrix(stimTypeArr == 1,:)).*100 ;
L_minus_M_Beta = mean(betaMatrix(stimTypeArr == 2,:)).*100 ;
S_Beta =         mean(betaMatrix(stimTypeArr == 3,:)).*100 ;

% Compute Standard Error
LightFluxBetaSE =  ((std(betaMatrix(stimTypeArr == 1,:)))./sqrt(numberOfRuns)).*100 ;
L_minus_M_BetaSE = ((std(betaMatrix(stimTypeArr == 2,:)))./sqrt(numberOfRuns)).*100 ;
S_BetaSE =         ((std(betaMatrix(stimTypeArr == 3,:)))./sqrt(numberOfRuns)).*100 ;

% Average Time Series for Each Combination of Stimulus Type & Run order
LightFluxAvgTS_A =  mean(timeSeriesStore(stimTypeArr == 1 & runOrder == 'A',:)) ;
L_minus_M_AvgTS_A = mean(timeSeriesStore(stimTypeArr == 2 & runOrder == 'A',:)) ;
S_AvgTS_A =         mean(timeSeriesStore(stimTypeArr == 3 & runOrder == 'A',:)) ;

LightFluxAvgTS_B =  mean(timeSeriesStore(stimTypeArr == 1 & runOrder == 'B',:)) ;
L_minus_M_AvgTS_B = mean(timeSeriesStore(stimTypeArr == 2 & runOrder == 'B',:)) ;
S_AvgTS_B =         mean(timeSeriesStore(stimTypeArr == 3 & runOrder == 'B',:)) ;

% Standard Error of Time Series
LightFluxStdTS_A =  (std(timeSeriesStore(stimTypeArr == 1 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;
L_minus_M_StdTS_A = (std(timeSeriesStore(stimTypeArr == 2 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;
S_StdTS_A =         (std(timeSeriesStore(stimTypeArr == 3 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder) ;

LightFluxStdTS_B =  (std(timeSeriesStore(stimTypeArr == 1 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;
L_minus_M_StdTS_B = (std(timeSeriesStore(stimTypeArr == 2 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;
S_StdTS_B =         (std(timeSeriesStore(stimTypeArr == 3 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder) ;

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
[wftd1, fp1] = fitWatsonToTTF_errorGuided(stimHz,LightFluxBeta,LightFluxBetaSE,1); hold on
errorbar(stimHz,LightFluxBeta,LightFluxBetaSE,'ko'); set(gca,'FontSize',15);
set(gca,'Xtick',stimHz); title('Light flux');

[wftd2, fp2] = fitWatsonToTTF_errorGuided(stimHz,L_minus_M_Beta,L_minus_M_BetaSE,1); hold on
errorbar(stimHz,L_minus_M_Beta,L_minus_M_BetaSE,'ko');
set(gca,'FontSize',15); set(gca,'Xtick',stimHz); title('L - M');
% S
[wftd3, fp3] = fitWatsonToTTF_errorGuided(stimHz,S_Beta,S_BetaSE,1); hold on
errorbar(stimHz,S_Beta,S_BetaSE,'ko'); set(gca,'FontSize',15);
set(gca,'Xtick',stimHz); title('S');

% Errorbars
figure;
errorbar(T_R.*(1:lengthAttnHRF)-1,mean(hrfStore),std(hrfStore)./sqrt(size(hrfStore,1)),'LineWidth',2);
xlabel('Time/s'); ylabel('% signal change'); title('HRF obtained using FIR');
set(gca,'FontSize',15);

%% Time Series plots 
% Use Function for plotting Data:
% -- plotLinModelFits -- 

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

