function f = forwardModelOutput(startTimesSorted,stimValuesSorted,tsFileNames, ...
          TS_timeSamples,stimDuration,stepFunctionRes,cosRamp, ...
          t_convolve,BOLDHRF,cleanedData,paramVec)

% function [betaMatrix,reconstructedTSmatrix,startTimesSorted_A, ...
%           stimValuesSorted_A, startTimesSorted_B, stimValuesSorted_B, ...
%           actualStimulusValues] ...
%                                       = ...
%           forwardModel(startTimesSorted,stimValuesSorted,tsFileNames, ...
%           TS_timeSamples,stimDuration,stepFunctionRes,cosRamp, ...
%           t_convolve,BOLDHRF,cleanedData)
%
% implements forward model for BOLD fitting project
%
% startTimesSorted: vector of starting times
% stimValuesSorted: stimulus values corresponding to startTimesSorted
% tsFileNames     : names of time series files corresponding to above
% TS_timesamples  : time samples for time series
% stimDuration    : stimulus duration
% stepFunctionRes : resolution for constructing stimulus vector
% cosRamp         : duration of cosine ramp, in seconds
% t_convolve      : time samples of convolution vectors
% cleanedData     : clean time series data

param = paramVec2structZhou(paramVec);
param.afterResponseTiming = 10;
param.MRamplitude = 1;
param.rectify = 1;

% get unique stimulus values
actualStimulusValues = unique(stimValuesSorted(stimValuesSorted~=-1 & stimValuesSorted~=0));

% Stire Stimulus Order A & B
stimValuesSorted_A = [] ;
stimValuesSorted_B = [] ;

% for each run
for i = 1:size(startTimesSorted,1)
    
    % Stores A & B sequences-- For labeling Time Series by Stimulus Period
   if strfind(char(tsFileNames(i)),'_A_') & isempty(stimValuesSorted_A)
      startTimesSorted_A = startTimesSorted(i,startTimesSorted(i,:)~=-1);  
      stimValuesSorted_A = stimValuesSorted(i,stimValuesSorted(i,:)~=-1); 
   elseif strfind(char(tsFileNames(i)),'_B_') & isempty(stimValuesSorted_B)
      startTimesSorted_B = startTimesSorted(i,startTimesSorted(i,:)~=-1); 
      stimValuesSorted_B = stimValuesSorted(i,stimValuesSorted(i,:)~=-1); 
   else
      stimOrderMarker = [] ; 
   end
   
    % for each stimulus
   for j = 1:length(actualStimulusValues)
       % grab the starting times
       startTimesForGivenStimValue = startTimesSorted(i,stimValuesSorted(i,:)==actualStimulusValues(j));
       % create stimulus model for each one, and sum those models together
       % to get the stimulus model for each stimulus type
       singleStimNeuralModel = [];
       for k = 1:length(startTimesForGivenStimValue)
           singleStimModel = createStimVector(TS_timeSamples,startTimesForGivenStimValue(k), ...
                        stimDuration,stepFunctionRes,cosRamp);
           yNeural = createNeuralTemporalModel(TS_timeSamples, singleStimModel, 0, param);
           singleStimNeuralModel(k,:) = yNeural;
       end
       % there is a vector for each run and stimulus
       singleStimModelNeuralAllRuns(i,j,:) = sum(singleStimNeuralModel);
   end
end

sumSquaredError = [];

% LOOP OVER RUNS
for i = 1:size(singleStimModelNeuralAllRuns,1)   
    designMatrixPreOnes = [];
    % LOOP OVER STIMULI WITHIN EACH RUN
   for j = 1:size(singleStimModelNeuralAllRuns,2)
       regressor = createRegressor(squeeze(singleStimModelNeuralAllRuns(i,j,:))',TS_timeSamples,BOLDHRF,t_convolve);
       % create design matrix (ones to be added later)
       designMatrixPreOnes(:,j) = regressor-mean(regressor);
   end
   % add ones regressor
   designMatrix = [ones([size(designMatrixPreOnes,1) 1]) designMatrixPreOnes] ;
    % Obtain Beta Weights
   betaWeights = designMatrix\cleanedData(i,:)' ; 

   % Beta Weights Sans Weight for the first Regressor
   betaMatrix(i,:) = betaWeights(2:length(betaWeights)) ;

   % Reconstruct Time Series According to Model 
   reconstructedTS = sum(repmat(betaWeights',[size(designMatrix,1) 1]).*designMatrix,2) ;

   % Store all reconstructed Time Series
   reconstructedTSmatrix(i,:) = reconstructedTS' ;
   
   sumSquaredError(i) = sum((reconstructedTS' - cleanedData(i,:)).^2);
end

f = sum(sumSquaredError);

gribble = 1;