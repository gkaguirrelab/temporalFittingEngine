function forwardModel(startTimesSorted,stimValuesSorted,tsFileNames, ...
                      TS_timeSamples,stimDuration,stepFunctionRes,cosRamp, ...
                      modelDuration,modelResolution,t_convolve,BOLDHRF, ...
                      singleStimModelAllRuns,cleanedData)

% function forwardModel(startTimesSorted,stimValuesSorted,tsFileNames)
%
% implements forward model for BOLD fitting project



gribble = 1;