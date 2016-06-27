function [betaMatrix,reconstructedTSmatrix,startTimesSorted_A, ...
          stimValuesSorted_A, startTimesSorted_B, stimValuesSorted_B, ...
          actualStimulusValues] ...
                                      = ...
          forwardModel(startTimesSorted,stimValuesSorted,tsFileNames, ...
          TS_timeSamples,stimDuration,stepFunctionRes,cosRamp, ...
          t_convolve,BOLDHRF,cleanedData)

% function forwardModel(startTimesSorted,stimValuesSorted,tsFileNames)
%
% implements forward model for BOLD fitting project

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
       singleStimModel = [];
       for k = 1:length(startTimesForGivenStimValue)
           singleStimModel(k,:) = createStimVector(TS_timeSamples,startTimesForGivenStimValue(k), ...
                        stimDuration,stepFunctionRes,cosRamp);
       end
       % there is a vector for each run and stimulus
       singleStimModelAllRuns(i,j,:) = sum(singleStimModel);
   end
end

% LOOP OVER RUNS
for i = 1:size(singleStimModelAllRuns,1)   
    designMatrixPreOnes = [];
    % LOOP OVER STIMULI WITHIN EACH RUN
   for j = 1:size(singleStimModelAllRuns,2)
       regressor = createRegressor(squeeze(singleStimModelAllRuns(i,j,:))',TS_timeSamples,BOLDHRF,t_convolve);
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
   reconstructedTSmatrix(i,:) = reconstructedTS ;
end

gribble = 1;