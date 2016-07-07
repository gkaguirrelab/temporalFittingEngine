function [stimMatrix,paramLockMatrix,startTimesSorted_A,startTimesSorted_B, ...
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
stimRef = [];
stimMatrix = [];
% Matrix for locking parameters
paramLockMatrix = [];

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
   
   startTimesForRun = startTimesSorted(i,stimValuesSorted(i,:)>0);
   stimValuesForRun = stimValuesSorted(i,stimValuesSorted(i,:)>0);
   
    % for each stimulus
   for j = 1:length(stimValuesForRun)
           % create the stimulus model
           stimVec = createStimVector(TS_timeSamples,startTimesForRun(j), ...
                        stimDuration,stepFunctionRes,cosRamp);
           stimMatrix(i,j,:) = stimVec; 
           assignmentMatrixRow = zeros([1 length(actualStimulusValues)]);
           stimRef(i,j,:) = double(stimValuesForRun(j) == actualStimulusValues)';           
   end
   % initialize for just this one run
   paramLockMatrixSubForCons1 = [];
   paramLockMatrixSub = [];
   % for each stimulus value
   for j = 1:length(actualStimulusValues)
      % find the indices at which that value occurs
      stimToBeChained = find(stimValuesForRun == actualStimulusValues(j));
      % then for each index, daisy-chain the linear equalities
      for k = 1:length(stimToBeChained)-1
         constraintVec = zeros([1 length(stimValuesForRun)]);
         constraintVec(stimToBeChained(k)) = 1;
         constraintVec(stimToBeChained(k+1)) = -1;
         paramLockMatrixSubForCons1(size(paramLockMatrixSubForCons1,1)+1,:) = constraintVec;
         paramLockMatrixSub = [paramLockMatrixSubForCons1 zeros(size(paramLockMatrixSubForCons1)); ...
                               zeros(size(paramLockMatrixSubForCons1)) paramLockMatrixSubForCons1];
         % testb = kron(eye(3),testa)
      end
   end
   % this only does the constraints for one run, so store it
   paramLockMatrix(i,:,:) = paramLockMatrixSub;
end
      
gribble = 1;