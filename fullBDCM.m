% Script for fitting 7T BOLD data elicited by cone contrast stimuli at
% different temporal frequencies. Fits with a non-linear model. Plots fits
% and fit parameters

%% Make sure toolboxes we need for this are on the path
AddToMatlabPathDynamically('BCDMToolbox');

%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_gka1' ; 
% *** Subject Pool ***
%     'HERO_asb1'--NOTE THAT ASB1 HAS 2 EXTRA LIGHTFLUX B AND 2 EXTRA S A
%                  RUNS--WE USE THE FIRST 2 OF EACH
%     'HERO_gka1'


session = 'all' ;
% *** collected on two dates ***
%     '041416' ...
%     '041516' ...
        
% 1 -> use canonical HRF, 0 -> extract HRF using FIR
bCanonicalHRF = 0;

% Boolean: 1 -> go into debug mode--only fit light flux A
bDEBUG = 0;

% do we turn off parameter locking?
bFreeFloatParams = 1;
%% LOAD TIME SERIES AND GET STIMULUS (& ATTENTION) START TIMES

% load time series
[avgTS, avgTSprc, tsFileNames, stimTypeCode, runOrder] ...
= loadTimeSeriesData(subj_name,session);

% get all stimulus values and start times, as well as the attention task
% start times
[startTimesSorted, stimValuesSorted, attnStartTimes] = orderStartTimes(subj_name,session);

% Time Series sampling points
timebase = [1:336]-1;

%% DERIVE HRF

% how long we expect the HRF to be in seconds
lengthHRF = 15;

% derive HRF from data
[BOLDHRF, cleanedData, SEHRF]= deriveHRFwrapper(avgTSprc,attnStartTimes,lengthHRF,'Fourier');

% downsample HRF
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

% deal with minor inconsistency in data collection files
if bFreeFloatParams == 1
   startTimesSorted = 0:stimDuration:max(timebase); 
   startTimesSorted = repmat(startTimesSorted,[size(stimValuesSorted,1) 1]);
end

% create stimulus vector
[stimMatrix,stimValues,startTimesSorted_A,startTimesSorted_B, ...
stimValuesSorted_A,stimValuesSorted_B,uniqueStimulusValues] ...
= createStimMatrix(startTimesSorted,stimValuesSorted,tsFileNames, ...
timebase,stimDuration,stepFunctionRes,cosRamp,bFreeFloatParams);

if bFreeFloatParams == 1
    stimValuesSorted_A = [0 stimValuesSorted_A];    
end

%% PRELIMINARIES FOR FITTING

% store parameters--initialize matrices
storeFitParams = zeros([0 0 0]);
storeAll = zeros([0 0 0]);

reconstructedTSmat = [];
errorStore = [];

% if in debug mode, only fit light flux A
if bDEBUG == 1
   runsToFit = find(stimTypeCode == 1 & runOrder == 'A');
else
   runsToFit = 1:size(stimMatrix,1);  
end
%% FIT
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
   if bFreeFloatParams == 1
      paramLockMatrix = []; 
   end
   % fit the response
   [paramsFit,fval,~,fittedResponse] = tmri.fitResponse({timebase},{stimulus},{boldResponse}, ...
                          'HRF',theHRF,'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
   fprintf('Model parameter from fits:\n');
   tmri.print(paramsFit);
   % store all parameters for each run, regardless of daisy-chaining
   storeAll(size(storeAll,1)+1,:,:) = paramsFit.paramMainMatrix;     
   % store RMS measure
   errorStore(length(errorStore)+1) = fval;
   [~,meanParamValues] = tmri.plotParams(paramsFit,stimulus,'bFitZero',logical(bFreeFloatParams)); close;
   % do this for each run
   storeFitParams(size(storeFitParams,1)+1,:,:) = meanParamValues;
   reconstructedTSmat(size(reconstructedTSmat,1)+1,:) = cell2mat(fittedResponse);
   display(['run number: ' num2str(i)]);
end
%% PLOTTING
if bDEBUG == 1
    % getting statistics over runs
   [meanMatrix, SEmatrix, stimTypeTagMatrix, paramNamesTagMatrix] = ...
    paramStatistics(permute(storeFitParams,[1 3 2]),stimTypeCode(runsToFit),paramsFit.paramNameCell);
    stimNamesCell = {'Light Flux'};    
    plotParamsWrapper(uniqueStimulusValues,meanMatrix,SEmatrix,stimTypeTagMatrix,paramNamesTagMatrix,stimNamesCell,[1 4 7])
    AvgTS = mean(cleanedData(runsToFit,:));
    StdTS =  std(cleanedData(runsToFit,:))./sqrt(size(storeFitParams,1));
    RMS = mean(errorStore);
    AvgTS_model = mean(reconstructedTSmat);
   
    % create cell for plotting stimulus starts
    stimValuesMatSorted_A_cell = {} ;
    for j = 1:length(stimValuesSorted_A)
       stimValuesMatSorted_A_cell{j} = num2str(stimValuesSorted_A(j)) ; 
    end
    % plot full time series
    figure;
    plotModelFits(timebase,AvgTS,AvgTS_model, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,StdTS,RMS);
    title('Light flux A'); xlabel('Time / s'); ylabel('% signal change');
elseif bDEBUG == 0
    % Parameter averaging    
    [meanMatrix, SEmatrix, stimTypeTagMatrix, paramNamesTagMatrix] = ...
    paramStatistics(permute(storeFitParams,[1 3 2]),stimTypeCode,paramsFit.paramNameCell);
    % TTF & tau2 plots
    stimNamesCell = {'Light Flux','L - M','S'};    
    plotParamsWrapper(uniqueStimulusValues,meanMatrix,SEmatrix,stimTypeTagMatrix,paramNamesTagMatrix,stimNamesCell)
    % get time series statistics
    [avgTS, stdTS, RMS, modelTS, idCell] = timeSeriesStatisticsBDCM(cleanedData,errorStore,reconstructedTSmat,stimTypeCode,runOrder);
    % plot
    plotModelFitsWrap(timebase, avgTS, stdTS, RMS, modelTS, idCell, ...
                         startTimesSorted_A, stimValuesSorted_A, ...
                         startTimesSorted_B, stimValuesSorted_B)
else
   error('fullBDCM: PICK VALID BOOLEAN VALUE FOR bDEBUG'); 
end

%% store everything in a struct in case we want to keep these plots
fitResults = struct; fitResults.storeFitParams = storeFitParams;
fitResults.stimTypeCode = stimTypeCode; fitResults.paramsFit = paramsFit;
fitResults.uniqueStimulusValues = uniqueStimulusValues; 
fitResults.cleanedData = cleanedData; fitResults.errorStore = errorStore;
fitResults.reconstructedTSmat = reconstructedTSmat; fitResults.runOrder = runOrder;
fitResults.timebase = timebase; fitResults.startTimesSorted_A = startTimesSorted_A;
fitResults.startTimesSorted_B = startTimesSorted_B;
fitResults.stimValuesSorted_A = stimValuesSorted_A;
fitResults.stimValuesSorted_B = stimValuesSorted_B;
fitResults.lengthHRF = lengthHRF; fitResults.hrf = hrf; 
fitResults.SEHRF = SEHRF; fitResults.hrfPointsToSample = hrfPointsToSample; 
fitResults.stimValues = stimValues; fitResults.stimValuesSorted = stimValuesSorted;
fitResults.startTimesSorted = startTimesSorted;
fitResults.subjName = subj_name; fitResults.storeAll = storeAll;
%% CARRY OVER MATRICES

% loop over runs
for i = 1:size(storeFitParams,1)
    stimForRun = stimValues(i,:);
    % get unique stimulus values
    uniqueTempFreq = unique(stimForRun);
    % initialize carry over matrices
    carryOverMatSubAmpl = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    carryOverMatSubtau2 = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    counterMatrixSub = zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    % loop over temp freqs in each run
    for j = 2:length(stimForRun)
        % get the amplitude and tau values
        ampl = storeAll(i,j,1);
        tau2 = storeAll(i,j,2);
        whereToPut1 = stimForRun(j)  ==uniqueTempFreq;
        whereToPut2 = stimForRun(j-1)==uniqueTempFreq;
        if isnan(carryOverMatSubAmpl(whereToPut1,whereToPut2))
            % place them appropriately in the carry over matrix for ampllitude
            carryOverMatSubAmpl(whereToPut1,whereToPut2) = ampl;
            % do the same for tau2
            carryOverMatSubtau2(whereToPut1,whereToPut2) = tau2;
        else
            % place them appropriately in the carry over matrix for amplitude
            carryOverMatSubAmpl(whereToPut1,whereToPut2) = ...
            carryOverMatSubAmpl(whereToPut1,whereToPut2)+ampl;
            % do the same for tau2
            carryOverMatSubtau2(whereToPut1,whereToPut2) = ...
            carryOverMatSubtau2(whereToPut1,whereToPut2)+tau2;
        end
        counterMatrixSub(whereToPut1,whereToPut2) = ...
        counterMatrixSub(whereToPut1,whereToPut2)+1;
    end
   
   % manually account for overlapping combinations
   if runOrder(i)=='A'
       
       carryOverMatSubAmpl(1,1) = carryOverMatSubAmpl(1,1)./2;
       carryOverMatSubAmpl(2,1) = carryOverMatSubAmpl(2,1)./2;
       
       carryOverMatSubtau2(1,1) = carryOverMatSubtau2(1,1)./2;
       carryOverMatSubtau2(2,1) = carryOverMatSubtau2(2,1)./2;
       
   elseif runOrder(i)=='B'
       
       carryOverMatSubAmpl(1,1) = carryOverMatSubAmpl(1,1)./3;
       carryOverMatSubtau2(1,1) = carryOverMatSubtau2(1,1)./3;
       
   else
      dummy = [];
   end
   
   % add 'labels'
   carryOverMatSubAmpl = [uniqueTempFreq' carryOverMatSubAmpl];
   carryOverMatSubAmpl = [[0 uniqueTempFreq]; carryOverMatSubAmpl];
   % store for each run
   carryOverMatAmpl(i,:,:) = carryOverMatSubAmpl;
   
   % repeat for tau2
   carryOverMatSubtau2 = [uniqueTempFreq' carryOverMatSubtau2];
   carryOverMatSubtau2 = [[0 uniqueTempFreq]; carryOverMatSubtau2];
   % store for each run
   carryOverMattau2(i,:,:) = carryOverMatSubtau2;
   
   % matrix for keeping track of repeats
   counterMatrix(i,:,:) = counterMatrixSub;
   
end

lightModDir = {'Light Flux','L - M','S'};

for i = 1:length(lightModDir)
   % create empty matrix
   finalCarryOverMatAmplSub = zeros(size(carryOverMatSubAmpl));
   finalCarryOverMattau2Sub = zeros(size(carryOverMatSubtau2));
   % get all runs with a given modulation direction
   runIndices = find(stimTypeCode==i);
   % grab their carry over matrices, sum them
   for j = 1:length(runIndices)
      % convert nan's to 0's
      curCarryOverAmpl = squeeze(carryOverMatAmpl(runIndices(j),:,:));
      curCarryOverTau2 = squeeze(carryOverMattau2(runIndices(j),:,:));
      curCarryOverAmpl(isnan(curCarryOverAmpl)) = 0;
      curCarryOverTau2(isnan(curCarryOverTau2)) = 0;
      % add
      finalCarryOverMatAmplSub = finalCarryOverMatAmplSub+curCarryOverAmpl; 
      finalCarryOverMattau2Sub = finalCarryOverMattau2Sub+curCarryOverTau2; 
   end
   % divide by half the number of runs at that modulation direction--two
   % run orders to get proper counterbalancing
   finalCarryOverMatAmplSub = finalCarryOverMatAmplSub./(length(runIndices)./2);
   finalCarryOverMattau2Sub = finalCarryOverMattau2Sub./(length(runIndices)./2);
   % divide the 'labels' by 6
   finalCarryOverMatAmplSub(1,:) = finalCarryOverMatAmplSub(1,:)./2;
   finalCarryOverMatAmplSub(:,1) = finalCarryOverMatAmplSub(:,1)./2;  
   finalCarryOverMattau2Sub(1,:) = finalCarryOverMattau2Sub(1,:)./2;
   finalCarryOverMattau2Sub(:,1) = finalCarryOverMattau2Sub(:,1)./2;
   % manually average overlapping cells
   finalCarryOverMatAmplSub(2,2) = finalCarryOverMatAmplSub(2,2)./2;
   finalCarryOverMattau2Sub(2,2) = finalCarryOverMattau2Sub(2,2)./2;
   finalCarryOverMatAmplSub(2,6) = finalCarryOverMatAmplSub(2,6)./2;
   finalCarryOverMattau2Sub(2,6) = finalCarryOverMattau2Sub(2,6)./2;
   % assign
   finalCarryOverMatAmpl(i,:,:) = finalCarryOverMatAmplSub;
   finalCarryOverMattau2(i,:,:) = finalCarryOverMattau2Sub;
   
end

for i = 1:length(lightModDir)
   plotCarryOver(uniqueTempFreq,squeeze(finalCarryOverMatAmpl( ...
                i,2:size(finalCarryOverMatAmpl,2),2:size(finalCarryOverMatAmpl,3))), ...
                'Amplitude',lightModDir(i))
   plotCarryOver(uniqueTempFreq,squeeze(finalCarryOverMattau2( ...
                i,2:size(finalCarryOverMattau2,2),2:size(finalCarryOverMattau2,3))), ...
                'tau2',lightModDir(i))
end

%% TRIAL AVERAGES

windowLength = 24;

% loop over stimulus modulation directions
for i = 1:length(lightModDir)
   % get all the unique temporal frequencies
   uniqueTempFreq = unique(stimValues(1,:));
   % find all indices for a given modulation direction
   runIndices = find(stimTypeCode==i);
   % get all data for given modulation direction
   cleanedDataForModDir = cleanedData(runIndices,:);
   % get corresponding fits
   reconTsForModDir = reconstructedTSmat(runIndices,:);
   % for each unique temporal frequency
   for j = 1:length(uniqueTempFreq)
       % initialize arrays for storing data and fits
       windowStoreData = [];
       windowStoreFit = [];
       % want to search across runs
       for k = 1:length(runIndices)
           % get the current run's stimulus values
           curRunOrder = stimValues(runIndices(k),:);
           % for each unique temporal frequency, get all times at which a
           % stimulus with its value came on
           startPoints = startTimesSorted(runIndices(k),uniqueTempFreq(j)==curRunOrder);
           % for each of these starting times
           for l = 1:length(startPoints)
              % grab all time points within the window, get their indices
              aWindow = find(timebase>=startPoints(l) ...
                        & timebase<=(startPoints(l)+windowLength));
              % store the data and reconstructed time series
              windowStoreData(size(windowStoreData,1)+1,1:length(aWindow)) = cleanedDataForModDir(k,aWindow);
              windowStoreFit(size(windowStoreFit,1)+1,1:length(aWindow)) = reconTsForModDir(k,aWindow);
           end
       end
       meanData = mean(windowStoreData);
       meanFit = mean(windowStoreFit);
       SEdata = std(windowStoreData)./sqrt(size(windowStoreData,1));
       SEfit = std(windowStoreFit)./sqrt(size(windowStoreFit,1));
       meanDataStore(i,j,:) = meanData;
       meanFitStore(i,j,:) = meanFit;
       SEdataStore(i,j,:) = SEdata;
       SEfitStore(i,j,:) = SEfit;
       display(num2str(size(windowStoreData,1)));
   end
end

figure;
set(gcf,'Position',[25 303 1872 788]);
for i = 1:size(meanDataStore,1)
    for j = 1:size(meanDataStore,2)
       subplot(length(lightModDir),length(uniqueTempFreq),(i-1).*length(uniqueTempFreq)+j) 
       plot(0:size(meanDataStore,3)-1,squeeze(meanDataStore(i,j,:))'); hold on
       plot(0:size(meanDataStore,3)-1,squeeze(meanFitStore(i,j,:))'); axis square;
       fill([0:size(meanDataStore,3)-1 fliplr(0:size(meanDataStore,3)-1)], ...
      [squeeze(meanDataStore(i,j,:))'+squeeze(SEdataStore(i,j,:))' ...
      fliplr(squeeze(meanDataStore(i,j,:))'-squeeze(SEdataStore(i,j,:))')], ...
      'k','FaceAlpha',0.15,'EdgeColor','none');
       if j == 1
          ylabel('% signal change');
       end
       if i == 3
          xlabel('time/s'); 
       end
       if j == 4 
          title([lightModDir(i) num2str(uniqueTempFreq(j)) ' Hz']);
       else
          title([num2str(uniqueTempFreq(j)) ' Hz']);
       end
    end
end