%%

subj_name = 'asb1';

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
    carryOverMatSubAmp = zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    carryOverMatSubtau2 = zeros([length(uniqueTempFreq) length(uniqueTempFreq)]);
    % loop over temp freqs in each run
    for j = 2:length(stimForRun)
        % get the amplitude and tau values
        amp = carryOver.storeAll(i,j,1);
        tau2 = carryOver.storeAll(i,j,2);
        % place them appropriately in the carry over matrix for amplitude
        carryOverMatSubAmp(stimForRun(j)==uniqueTempFreq, ...
        stimForRun(j-1)==uniqueTempFreq) = amp;
        % do the same for tau2
        carryOverMatSubtau2(stimForRun(j)==uniqueTempFreq, ...
        stimForRun(j-1)==uniqueTempFreq) = tau2;
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
      finalCarryOverMatAmpSub = finalCarryOverMatAmpSub+squeeze(carryOverMatAmp(runIndices(j),:,:)); 
      finalCarryOverMattau2Sub = finalCarryOverMattau2Sub+squeeze(carryOverMattau2(runIndices(j),:,:)); 
   end
   % divide by half the number of runs at that modulation direction--two
   % run orders to get proper counterbalancing
   finalCarryOverMatAmpSub = finalCarryOverMatAmpSub./(length(runIndices)./2);
   finalCarryOverMattau2Sub = finalCarryOverMattau2Sub./(length(runIndices)./2);
   % divide the 'labels' by 2
   finalCarryOverMatAmpSub(1,:) = finalCarryOverMatAmpSub(1,:)./2;
   finalCarryOverMatAmpSub(:,1) = finalCarryOverMatAmpSub(:,1)./2;
   finalCarryOverMattau2Sub(1,:) = finalCarryOverMattau2Sub(1,:)./2;
   finalCarryOverMattau2Sub(:,1) = finalCarryOverMattau2Sub(:,1)./2;
   % assign
   finalCarryOverMatAmp(i,:,:) = finalCarryOverMatAmpSub;
   finalCarryOverMattau2(i,:,:) = finalCarryOverMattau2Sub;
end