function [avgTS, stdTS, MSE, modelTS, idCell] = paramStatisticsBDCM(cleanedData,MSEstore,reconstructedTSmat,stimTypeArr,runOrder)

uniqueStimType = unique(stimTypeArr);
uniqueRunOrder = unique(runOrder);

avgTS = [];
stdTS = [];
MSE = [];
modelTS = [];
modDirCell = {'Light Flux','L - M','S'};
idCell = {};

for i = 1:length(uniqueStimType)
    for j = 1:length(uniqueRunOrder)
       curInd = find(stimTypeArr == uniqueStimType(i) & runOrder == uniqueRunOrder(j));
       curCleanedData = cleanedData(curInd,:);
       avgTS(size(avgTS,1)+1,:) = mean(curCleanedData);
       stdTS(size(stdTS,1)+1,:) = std(curCleanedData);
       MSE(size(MSE,1)+1,:) = mean(MSEstore(curInd));
       curReconstructedTS = reconstructedTSmat(curInd,:);
       modelTS(size(modelTS,1)+1,:) = mean(curReconstructedTS);
       idRun = [char(modDirCell(uniqueStimType(i))) ' ' uniqueRunOrder(j)];
       idCell{length(idCell)+1} = idRun;
    end
end

end