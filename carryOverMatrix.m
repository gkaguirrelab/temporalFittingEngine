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
    stimForRun = carryOver.stimValuesSorted(i,carryOver.stimValuesSorted(i,:)~=-1);
    if carryOver.runOrder(i) == 'A'
       stimForRun = [0 stimForRun]; 
    end
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
        carryOverMatSubAmp(carryOver.stimValuesForRunStore(i,j)==uniqueTempFreq, ...
        carryOver.stimValuesForRunStore(i,j-1)==uniqueTempFreq) = amp;
        % do the same for tau2
        carryOverMatSubtau2(carryOver.stimValuesForRunStore(i,j)==uniqueTempFreq, ...
        carryOver.stimValuesForRunStore(i,j-1)==uniqueTempFreq) = tau2;
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