%%

subj_name = 'gka1';
% NOTE THAT THE COUNTERBALANCING ISN'T PERFECT--MANUALLY CORRECTED IN CODE
% FOR RUN A: 2 PRECEDED BY 0 APPEARS 2 TIMES, 0 PRECEDED BY 0 APPEARS 2
%            TIMES
% FOR RUN B: 0 PRECEDED BY 0 APPEARS 3 TIMES
hardCodedOverlapCorrection = 2;

% load parameters and fits
if strcmp(subj_name,'asb1')
    carryOver = load('asb1CarryOver');
    carryOver = carryOver.carryOver;
    unlockedParamFits = load('asb1UnlockedParamFits');
    unlockedParamFits = unlockedParamFits.unlockedParamFits;
elseif strcmp(subj_name,'gka1')
    carryOver = load('gka1CarryOver');
    carryOver = carryOver.carryOver;
    unlockedParamFits = load('gka1UnlockedParamFits');
    unlockedParamFits = unlockedParamFits.unlockedParamFits;
else
   error('specify valid subject name: asb1 or gka1'); 
end

% loop over runs
for i = 1:size(carryOver.storeAll,1)
    stimForRun = carryOver.stimValuesForRunStore(i,:);
    % get unique stimulus values
    uniqueTempFreq = unique(stimForRun);
    % initialize carry over matrices
    carryOverMatSubAmp = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    carryOverMatSubtau2 = 0./zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    counterMatrixSub = zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    % loop over temp freqs in each run
    for j = 2:length(stimForRun)
        % get the amplitude and tau values
        amp = carryOver.storeAll(i,j,1);
        tau2 = carryOver.storeAll(i,j,2);
        whereToPut1 = stimForRun(j)==uniqueTempFreq;
        whereToPut2 = stimForRun(j-1)==uniqueTempFreq;
        if isnan(carryOverMatSubAmp(whereToPut1,whereToPut2))
            % place them appropriately in the carry over matrix for amplitude
            carryOverMatSubAmp(whereToPut1,whereToPut2) = amp;
            % do the same for tau2
            carryOverMatSubtau2(whereToPut1,whereToPut2) = tau2;
        else
            % place them appropriately in the carry over matrix for amplitude
            carryOverMatSubAmp(whereToPut1,whereToPut2) = ...
            carryOverMatSubAmp(whereToPut1,whereToPut2)+amp;
            % do the same for tau2
            carryOverMatSubtau2(whereToPut1,whereToPut2) = ...
            carryOverMatSubAmp(whereToPut1,whereToPut2)+tau2;
        end
        counterMatrixSub(whereToPut1,whereToPut2) = ...
        counterMatrixSub(whereToPut1,whereToPut2)+1;
    end
   
   % manually account for overlapping combinations
   if carryOver.runOrder(i)=='A'
       carryOverMatSubAmp(2,2) = carryOverMatSubAmp(2,2)./2;
       carryOverMatSubAmp(3,2) = carryOverMatSubAmp(3,2)./2;
       carryOverMatSubtau2(2,2) = carryOverMatSubtau2(2,2)./2;
       carryOverMatSubtau2(3,2) = carryOverMatSubtau2(3,2)./2;
   elseif carryOver.runOrder(i)=='B'
       carryOverMatSubAmp(2,2) = carryOverMatSubAmp(2,2)./3;
       carryOverMatSubtau2(2,2) = carryOverMatSubtau2(2,2)./3;
   else
      dummy = [];
   end
   % add 'labels'
   carryOverMatSubAmp = [uniqueTempFreq' carryOverMatSubAmp];
   carryOverMatSubAmp = [[0 uniqueTempFreq]; carryOverMatSubAmp];
   % store for each run
   carryOverMatAmp(i,:,:) = carryOverMatSubAmp;
   
   % repeat for tau2
   carryOverMatSubtau2 = [uniqueTempFreq' carryOverMatSubtau2];
   carryOverMatSubtau2 = [[0 uniqueTempFreq]; carryOverMatSubtau2];
   carryOverMattau2(i,:,:) = carryOverMatSubtau2;
   % matrix for keeping track of repeats
   counterMatrix(i,:,:) = counterMatrixSub;
   
end

lightModDir = {'Light Flux','L - M','S'};

for i = 1:length(lightModDir)
   % create empty matrix
   finalCarryOverMatAmpSub = zeros(size(carryOverMatSubAmp));
   finalCarryOverMattau2Sub = zeros(size(carryOverMatSubtau2));
   % get all runs with a given modulation direction
   runIndices = find(carryOver.stimTypeCode==i);
   % grab their carry over matrices, sum them
   for j = 1:length(runIndices)
      % convert nan's to 0's
      curCarryOverAmp = squeeze(carryOverMatAmp(runIndices(j),:,:));
      curCarryOverTau2 = squeeze(carryOverMattau2(runIndices(j),:,:));
      curCarryOverAmp(isnan(curCarryOverAmp)) = 0;
      curCarryOverTau2(isnan(curCarryOverTau2)) = 0;
      % add
      finalCarryOverMatAmpSub = finalCarryOverMatAmpSub+curCarryOverAmp; 
      finalCarryOverMattau2Sub = finalCarryOverMattau2Sub+curCarryOverTau2; 
   end
   % divide by half the number of runs at that modulation direction--two
   % run orders to get proper counterbalancing
   finalCarryOverMatAmpSub = finalCarryOverMatAmpSub./(length(runIndices)./6);
   finalCarryOverMattau2Sub = finalCarryOverMattau2Sub./(length(runIndices)./6);
   % divide the 'labels' by 2
   finalCarryOverMatAmpSub(1,:) = finalCarryOverMatAmpSub(1,:)./2;
   finalCarryOverMatAmpSub(:,1) = finalCarryOverMatAmpSub(:,1)./2;
   finalCarryOverMattau2Sub(1,:) = finalCarryOverMattau2Sub(1,:)./2;
   finalCarryOverMattau2Sub(:,1) = finalCarryOverMattau2Sub(:,1)./2;
   % manually average overlapping cells
   finalCarryOverMatAmpSub(2,2) = finalCarryOverMatAmpSub(2,2)./2;
   finalCarryOverMattau2Sub(2,2) = finalCarryOverMattau2Sub(2,2)./2;
   finalCarryOverMatAmpSub(2,6) = finalCarryOverMatAmpSub(2,6)./2;
   finalCarryOverMattau2Sub(2,6) = finalCarryOverMattau2Sub(2,6)./2;
   % assign
   finalCarryOverMatAmp(i,:,:) = finalCarryOverMatAmpSub;
   finalCarryOverMattau2(i,:,:) = finalCarryOverMattau2Sub;
end


figure; imagesc(squeeze(finalCarryOverMatAmp(1,2:size(finalCarryOverMatAmp,2),2:size(finalCarryOverMatAmp,3))))
colormap gray; title('Light Flux'); set(gca,'xticklabel',([0 2 4 8 16 32 64]));
set(gca,'yticklabel',([0 2 4 8 16 32 64])); xlabel('Preceding stimulus (Hz)');
ylabel('Stimulus'); set(gca,'FontSize',15);

figure; imagesc(squeeze(finalCarryOverMatAmp(2,2:size(finalCarryOverMatAmp,2),2:size(finalCarryOverMatAmp,3))))
colormap gray; title('L - M'); set(gca,'xticklabel',([0 2 4 8 16 32 64]));
set(gca,'yticklabel',([0 2 4 8 16 32 64])); xlabel('Preceding stimulus (Hz)');
ylabel('Stimulus'); set(gca,'FontSize',15);

figure; imagesc(squeeze(finalCarryOverMatAmp(3,2:size(finalCarryOverMatAmp,2),2:size(finalCarryOverMatAmp,3))))
colormap gray; title('S'); set(gca,'xticklabel',([0 2 4 8 16 32 64]));
set(gca,'yticklabel',([0 2 4 8 16 32 64])); xlabel('Preceding stimulus (Hz)');
ylabel('Stimulus'); set(gca,'FontSize',15);