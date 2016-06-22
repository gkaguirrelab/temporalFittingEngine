function plotLinModelFits(timeSamples,timeSeriesAvgAct,timeSeriesAvgModel,startTimes,stimValuesCell,stimValuesMat,timeSeriesStd)

% function plotLinModelFits(timeSamples,timeSeriesAvgAct,timeSeriesAvgModel,startTimes,stimValuesCell,stimValuesMat,timeSeriesStd)
%
% specialized function for making plots

plot(timeSamples,timeSeriesAvgAct); hold on
plot(timeSamples,timeSeriesAvgModel);
text(startTimes,repmat(max(timeSeriesAvgAct),[1 length(startTimes)]),stimValuesCell);
fill([timeSamples fliplr(timeSamples)], ...
     [timeSeriesAvgAct+timeSeriesStd fliplr(timeSeriesAvgAct-timeSeriesStd)],'k','FaceAlpha',0.15,'EdgeColor','none');
makeStimColorLine(startTimes, ...
                  repmat(min(timeSeriesAvgAct-timeSeriesStd), ...
                  [1 length(startTimes)]),stimValuesMat)

gribble = 1;