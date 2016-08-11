%% MAKES FULL SET OF BDCM PLOTS

% make sure to add mriTemporalFitting and all subdirectories to path. that
% is, cd into mriTemporalFitting, then call addpath(genpath(pwd))

subjName = 'gka1';

if strcmp(subjName,'gka1')
    load('gka1_fitResults');
elseif strcmp(subjName,'asb1')
    load('asb1_fitResults');
else
   error('makeFullSetOfPlots: pick valid subject'); 
end

figure;
errorbar(0:fitResults.lengthHRF,fitResults.hrf,fitResults.SEHRF(fitResults.hrfPointsToSample),'LineWidth',2)
xlabel('Time/s'); ylabel('Signal'); set(gca,'FontSize',15);
title('HRF');

% Parameter averaging    
[meanMatrix, SEmatrix, stimTypeTagMatrix, paramNamesTagMatrix] = ...
paramStatistics(permute(fitResults.storeFitParams,[1 3 2]),fitResults.stimTypeCode,fitResults.paramsFit.paramNameCell);
% TTF & tau2 plots
stimNamesCell = {'Light Flux','L - M','S'};    
plotParamsWrapper(fitResults.uniqueStimulusValues,meanMatrix,SEmatrix,stimTypeTagMatrix,paramNamesTagMatrix,stimNamesCell)
% get time series statistics
[avgTS, stdTS, RMS, modelTS, idCell] = timeSeriesStatisticsBDCM(fitResults.cleanedData,fitResults.errorStore,fitResults.reconstructedTSmat,fitResults.stimTypeCode,fitResults.runOrder);
% plot
plotModelFitsWrap(fitResults.timebase, avgTS, stdTS, RMS, modelTS, idCell, ...
                     fitResults.startTimesSorted_A, fitResults.stimValuesSorted_A, ...
                     fitResults.startTimesSorted_B, fitResults.stimValuesSorted_B)
% 
% % loop over runs
% for i = 1:size(fitResults.storeFitParams,1)
%     stimForRun = fitResults.stimValues(i,:);
%     % get unique stimulus values
%     uniqueTempFreq = unique(stimForRun);
%     % initialize carry over matrices
%     carryOverMatSubAmpl = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
%     carryOverMatSubtau2 = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
%     counterMatrixSub = zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
%     % loop over temp freqs in each run
%     for j = 2:length(stimForRun)
%         % get the amplitude and tau values
%         ampl = fitResults.storeFitParams(i,1,j);
%         tau2 = fitResults.storeFitParams(i,2,j);
%         whereToPut1 = stimForRun(j)  ==uniqueTempFreq;
%         whereToPut2 = stimForRun(j-1)==uniqueTempFreq;
%         if isnan(carryOverMatSubAmpl(whereToPut1,whereToPut2))
%             % place them appropriately in the carry over matrix for ampllitude
%             carryOverMatSubAmpl(whereToPut1,whereToPut2) = ampl;
%             % do the same for tau2
%             carryOverMatSubtau2(whereToPut1,whereToPut2) = tau2;
%         else
%             % place them appropriately in the carry over matrix for amplitude
%             carryOverMatSubAmpl(whereToPut1,whereToPut2) = ...
%             carryOverMatSubAmpl(whereToPut1,whereToPut2)+ampl;
%             % do the same for tau2
%             carryOverMatSubtau2(whereToPut1,whereToPut2) = ...
%             carryOverMatSubtau2(whereToPut1,whereToPut2)+tau2;
%         end
%         counterMatrixSub(whereToPut1,whereToPut2) = ...
%         counterMatrixSub(whereToPut1,whereToPut2)+1;
%     end
%    
%    % manually account for overlapping combinations
%    if fitResults.runOrder(i)=='A'
%        
%        carryOverMatSubAmpl(1,1) = carryOverMatSubAmpl(1,1)./2;
%        carryOverMatSubAmpl(2,1) = carryOverMatSubAmpl(2,1)./2;
%        
%        carryOverMatSubtau2(1,1) = carryOverMatSubtau2(1,1)./2;
%        carryOverMatSubtau2(2,1) = carryOverMatSubtau2(2,1)./2;
%        
%    elseif fitResults.runOrder(i)=='B'
%        
%        carryOverMatSubAmpl(1,1) = carryOverMatSubAmpl(1,1)./3;
%        carryOverMatSubtau2(1,1) = carryOverMatSubtau2(1,1)./3;
%        
%    else
%       dummy = [];
%    end
%    
%    % add 'labels'
%    carryOverMatSubAmpl = [uniqueTempFreq' carryOverMatSubAmpl];
%    carryOverMatSubAmpl = [[0 uniqueTempFreq]; carryOverMatSubAmpl];
%    % store for each run
%    carryOverMatAmpl(i,:,:) = carryOverMatSubAmpl;
%    
%    % repeat for tau2
%    carryOverMatSubtau2 = [uniqueTempFreq' carryOverMatSubtau2];
%    carryOverMatSubtau2 = [[0 uniqueTempFreq]; carryOverMatSubtau2];
%    % store for each run
%    carryOverMattau2(i,:,:) = carryOverMatSubtau2;
%    
%    % matrix for keeping track of repeats
%    counterMatrix(i,:,:) = counterMatrixSub;
%    
% end
% 
% lightModDir = {'Light Flux','L - M','S'};
% 
% for i = 1:length(lightModDir)
%    % create empty matrix
%    finalCarryOverMatAmplSub = zeros(size(carryOverMatSubAmpl));
%    finalCarryOverMattau2Sub = zeros(size(carryOverMatSubtau2));
%    % get all runs with a given modulation direction
%    runIndices = find(fitResults.stimTypeCode==i);
%    % grab their carry over matrices, sum them
%    for j = 1:length(runIndices)
%       % convert nan's to 0's
%       curCarryOverAmpl = squeeze(carryOverMatAmpl(runIndices(j),:,:));
%       curCarryOverTau2 = squeeze(carryOverMattau2(runIndices(j),:,:));
%       curCarryOverAmpl(isnan(curCarryOverAmpl)) = 0;
%       curCarryOverTau2(isnan(curCarryOverTau2)) = 0;
%       % add
%       finalCarryOverMatAmplSub = finalCarryOverMatAmplSub+curCarryOverAmpl; 
%       finalCarryOverMattau2Sub = finalCarryOverMattau2Sub+curCarryOverTau2; 
%    end
%    % divide by half the number of runs at that modulation direction--two
%    % run orders to get proper counterbalancing
%    finalCarryOverMatAmplSub = finalCarryOverMatAmplSub./(length(runIndices)./2);
%    finalCarryOverMattau2Sub = finalCarryOverMattau2Sub./(length(runIndices)./2);
%    % divide the 'labels' by 6
%    finalCarryOverMatAmplSub(1,:) = finalCarryOverMatAmplSub(1,:)./2;
%    finalCarryOverMatAmplSub(:,1) = finalCarryOverMatAmplSub(:,1)./2;
%    finalCarryOverMattau2Sub(1,:) = finalCarryOverMattau2Sub(1,:)./2;
%    finalCarryOverMattau2Sub(:,1) = finalCarryOverMattau2Sub(:,1)./2;
%    % manually average overlapping cells
%    finalCarryOverMatAmplSub(2,2) = finalCarryOverMatAmplSub(2,2)./2;
%    finalCarryOverMattau2Sub(2,2) = finalCarryOverMattau2Sub(2,2)./2;
%    finalCarryOverMatAmplSub(2,6) = finalCarryOverMatAmplSub(2,6)./2;
%    finalCarryOverMattau2Sub(2,6) = finalCarryOverMattau2Sub(2,6)./2;
%    % assign
%    finalCarryOverMatAmpl(i,:,:) = finalCarryOverMatAmplSub;
%    finalCarryOverMattau2(i,:,:) = finalCarryOverMattau2Sub;
%    
% end