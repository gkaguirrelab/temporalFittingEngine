%%

stimNamesCell = {'Light Flux','L - M','S'};

[meanMatrix, SEmatrix, stimTypeTagMatrix, paramNamesTagMatrix] = ...
paramStatistics(storeUnique,stimTypeArr,paramStructFit.paramNameCell);

plotParamsWrapper(actualStimulusValues,meanMatrix,SEmatrix,stimTypeTagMatrix,paramNamesTagMatrix,stimNamesCell)