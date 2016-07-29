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

lightModDir = {'Light Flux','L - M','S'};

T_R = 1;

timeBase = 0:T_R:T_R*length(unlockedParamFits.cleanedData(1,:))-1;

windowLength = 24;

% loop over stimulus modulation directions
for i = 1:length(lightModDir)
   % get all the unique temporal frequencies
   uniqueTempFreq = unique(carryOver.stimValuesForRunStore(1,:));
   % find all indices for a given modulation direction
   runIndices = find(carryOver.stimTypeCode==i);
   % get all data for given modulation direction
   cleanedDataForModDir = unlockedParamFits.cleanedData(runIndices,:);
   % get corresponding fits
   reconTsForModDir = unlockedParamFits.reconstructedTSmat(runIndices,:);
   % for each unique temporal frequency
   for j = 1:length(uniqueTempFreq)
       % initialize arrays for storing data and fits
       windowStoreData = [];
       windowStoreFit = [];
       % want to search across runs
       for k = 1:length(runIndices)
           % get the current run's stimulus values
           curRunOrder = carryOver.stimValuesForRunStore(runIndices(k),:);
           % for each unique temporal frequency, get all times at which a
           % stimulus with its value came on
           startPoints = carryOver.startTimesSorted(runIndices(k),uniqueTempFreq(j)==curRunOrder);
           % for each of these starting times
           for l = 1:length(startPoints)
              % grab all time points within the window, get their indices
              aWindow = find(timeBase>=startPoints(l) ...
                        & timeBase<=(startPoints(l)+windowLength));
              % store the data and reconstructed time series
              windowStoreData(size(windowStoreData,1)+1,1:length(aWindow)) = cleanedDataForModDir(k,aWindow);
              windowStoreFit(size(windowStoreFit,1)+1,1:length(aWindow)) = reconTsForModDir(k,aWindow);
           end
       end
   end
end
