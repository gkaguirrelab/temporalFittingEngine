function [stimMatrix,stimValuesForRunStore,startTimesSorted_A,startTimesSorted_B, ...
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
stimMatrix = [];
% Matrix for locking parameters
stimValuesForRunStore = [];

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
   
   % start times and corresponding stim values--get rid of filler numbers
   startTimesForRun = startTimesSorted(i,stimValuesSorted(i,:)>0);
   stimValuesForRun = stimValuesSorted(i,stimValuesSorted(i,:)>0);
   % store for param locking in main function
   stimValuesForRunStore(i,:) = stimValuesForRun;
   
    % for each stimulus
   for j = 1:length(stimValuesForRun)
           % create the stimulus model
           stimVec = createStimVector(TS_timeSamples,startTimesForRun(j), ...
                        stimDuration,stepFunctionRes,cosRamp);
           stimMatrix(i,j,:) = stimVec; 
   end
   
end
      
gribble = 1;