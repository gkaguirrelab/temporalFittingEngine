function [stimMatrix,startTimesSorted_A,startTimesSorted_B, ...
          stimValuesSorted_A,stimValuesSorted_B,actualStimulusValues] ...
          = createStimMatrix(startTimesSorted,stimValuesSorted, ...
          tsFileNames,TS_timeSamples,stimDuration,stepFunctionRes,cosRamp)
      
% function stimMatrix = createStimMatrix(startTimesSorted,stimValuesSorted, ...
%           TS_timeSamples,stimDuration,stepFunctionRes,cosRamp)
%
% creates stimulus matrix

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
           % create the stimulus model
           singleStimModel(k,:) = createStimVector(TS_timeSamples,startTimesForGivenStimValue(k), ...
                        stimDuration,stepFunctionRes,cosRamp);
       end
       stimMatrix(i,j,:) = sum(singleStimModel);       
   end
end
      
gribble = 1;