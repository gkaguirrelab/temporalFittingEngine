%% temporalTransferFunc

% CREATES TEMPORAL TRANSFER FUNCTIONS.  

% SPECIFY THE SUBJECT AND DATE BELOW. 

%% Variable name legend (will clean these up if there is time)

% ts = time series
% LH and RH = left and right hemisphere
% stim = stimulus
% AVG = average
% Arr = array (superfluous--one of Ben's coding habit)
% Mat = matrix
% cur = current
% position = index

%% Identify the user
 if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
 end


%% SPECIFY SUBJECT AND SESSION, AND DROPBOX FOLDER

subj_name = 'HERO_gka1';
%     'HERO_asb1' 
%     'HERO_gka1'


session = 'all';
%     '041416' ...
%     '041516' ...

% PATH TO LOCAL DROPBOX
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'];

%% HRF PARAMETERS (GRABBED FROM WINAWER MODEL CODE)

% 1 = USE CANONICAL HRF, 0 = USE THE FIR-EXTRACTED HRF
bCanonicalHRF = 0;

if bCanonicalHRF == 1
    param = struct;
    % parameters of the double-gamma hemodynamic filter (HRF)
    param.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in seconds)
    param.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in seconds)
    param.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
end

% DURATION OF STIMULUS (ALWAYS THE SAME)
stimTime = 12;
% DURATION OF ATTENTION TASK
attnTime = 0.25;
% THE ATTENTION TASK HAS NO FREQUENCY, BUT IT GETS A 'CODE'
attnCode = 1;
% BOOLEAN FOR WHETHER WE WANT TO MODEL THE ATTENTION TASK WITH FIR
attnFIR = 1;
% LENGTH OF ATTENTION HRF
lengthAttnHRF = 26;
% ACQUISITION TIME
T_R = 1;
% TIME SERIES SAMPLING POINTS
TS_timeSamples = 1:336;
        
%% DEFINING PATHS, ORDER, ETC.

% DEFINE PATH TO FOLDER FOR ONE SUBJECT ON ONE DATE
dirPathStim = [localDropboxDir 'MELA_analysis/HCLV_Photo_7T/mriTemporalFitting_data/' ...
           subj_name '/' session '/' 'Stimuli/'];

% DEFINE PATH TO TIME SERIES DATA
dirPathTimeSeries = [localDropboxDir 'MELA_analysis/HCLV_Photo_7T/mriTemporalFitting_data/' ...
                      subj_name '/' 'TimeSeries/'];

% ORDER IN WHICH TIME SERIES DATA WAS COLLECTED, FOR FIGURING OUT WHICH
% TIME SERIES' TO PLOT WITH WHICH STIMULUS
tsFileNamesASB1_DAY1 = { ...
    'bold_1.6_P2_mb5_LightFlux_A_run1' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run1' ...
    'bold_1.6_P2_mb5_S_A_run1' ...
    'bold_1.6_P2_mb5_S_B_run1' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run1' ...
    'bold_1.6_P2_mb5_LightFlux_A_run2' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run2' ...
    'bold_1.6_P2_mb5_S_B_run2' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run2' ...
    'bold_1.6_P2_mb5_LightFlux_B_run2'
};

tsFileNamesASB1_DAY2 = { ...    
    'bold_1.6_P2_mb5_LightFlux_A_run3' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run3' ...
    'bold_1.6_P2_mb5_S_A_run3' ...
    'bold_1.6_P2_mb5_S_B_run3' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run3' ...
    'bold_1.6_P2_mb5_LightFlux_B_run3' ...
    'bold_1.6_P2_mb5_LightFlux_A_run4' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run4' ...
    'bold_1.6_P2_mb5_S_A_run4' ...
    'bold_1.6_P2_mb5_S_B_run4' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run4' ...
    'bold_1.6_P2_mb5_LightFlux_B_run4' ...   
    'bold_1.6_P2_mb5_LightFlux_A_run5' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run5' ...
    'bold_1.6_P2_mb5_S_A_run5' ...    
    'bold_1.6_P2_mb5_S_B_run5' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run5' ...
    'bold_1.6_P2_mb5_LightFlux_B_run5' ...
    'bold_1.6_P2_mb5_LightFlux_A_run6' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run6' ...
    'bold_1.6_P2_mb5_S_A_run6' ...
    'bold_1.6_P2_mb5_S_B_run6' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run6' ...
    'bold_1.6_P2_mb5_LightFlux_B_run6' ...
    'bold_1.6_P2_mb5_LightFlux_B_run1' ...
    'bold_1.6_P2_mb5_S_A_run2' ...
};

tsFileNamesGKA1_DAY1 = { ...
    'bold_1.6_P2_mb5_LightFlux_A_run1' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run1' ...
    'bold_1.6_P2_mb5_S_A_run1' ...
    'bold_1.6_P2_mb5_S_B_run1' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run1' ...
    'bold_1.6_P2_mb5_LightFlux_B_run1' ...
    'bold_1.6_P2_mb5_LightFlux_A_run2' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run2' ...
    'bold_1.6_P2_mb5_S_A_run2' ...
    'bold_1.6_P2_mb5_S_B_run2' ...    
    'bold_1.6_P2_mb5_L_minus_M_B_run2' ...
    'bold_1.6_P2_mb5_LightFlux_B_run2' ...
    'bold_1.6_P2_mb5_LightFlux_A_run3' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run3' ...
    'bold_1.6_P2_mb5_S_A_run3' ...
    'bold_1.6_P2_mb5_S_B_run3' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run3' ...
    'bold_1.6_P2_mb5_LightFlux_B_run3' ...
    'bold_1.6_P2_mb5_LightFlux_A_run4' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run4' ...
    'bold_1.6_P2_mb5_S_A_run4' ...
};

tsFileNamesGKA1_DAY2 = { ...
    'bold_1.6_P2_mb5_S_B_run4' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run4' ...
    'bold_1.6_P2_mb5_LightFlux_B_run4' ...
    'bold_1.6_P2_mb5_LightFlux_A_run5' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run5' ...
    'bold_1.6_P2_mb5_S_A_run5' ...
    'bold_1.6_P2_mb5_S_B_run5' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run5' ...
    'bold_1.6_P2_mb5_LightFlux_B_run5' ...
    'bold_1.6_P2_mb5_LightFlux_A_run6' ...
    'bold_1.6_P2_mb5_L_minus_M_A_run6' ...
    'bold_1.6_P2_mb5_S_A_run6' ...
    'bold_1.6_P2_mb5_S_B_run6' ...
    'bold_1.6_P2_mb5_L_minus_M_B_run6' ...
    'bold_1.6_P2_mb5_LightFlux_B_run6' ...
};

%% GET TIME SERIES DATA

% SUBJECT AND DATE DETERMINE WHICH TIME SERIES FILES WE LOAD, AND THE ORDER
% THEY ARE PLOTTED IN
if strcmp(subj_name,'HERO_asb1') & strcmp(session,'041416')
    currentTimeSeriesFolder = tsFileNamesASB1_DAY1;
elseif strcmp(subj_name,'HERO_asb1') & strcmp(session,'041516')
    currentTimeSeriesFolder = tsFileNamesASB1_DAY2;
elseif strcmp(subj_name,'HERO_asb1') & strcmp(session,'all')
    currentTimeSeriesFolder = {tsFileNamesASB1_DAY1{:},tsFileNamesASB1_DAY2{:}};
elseif strcmp(subj_name,'HERO_gka1') & strcmp(session,'041416')
    currentTimeSeriesFolder = tsFileNamesGKA1_DAY1;
elseif strcmp(subj_name,'HERO_gka1') & strcmp(session,'041516')
    currentTimeSeriesFolder = tsFileNamesGKA1_DAY2;
elseif strcmp(subj_name,'HERO_gka1') & strcmp(session,'all')
    currentTimeSeriesFolder = {tsFileNamesGKA1_DAY1{:},tsFileNamesGKA1_DAY2{:}};
else
    error('stimulusTimeSeries ERROR: INPUT SUBJECT / DATE DOES NOT EXIST');
end

% GET THE CONTENTS OF THE TIME SERIES FOLDER
timeSeriesDir = dir(dirPathTimeSeries);

% GET A CELL CONTAINING THEIR NAMES BY LOOPING OVER
timeSeriesDirNames = {};

for i = 1:length(timeSeriesDir)
   tsFileName = timeSeriesDir(i).name;
   if length(tsFileName)>15
       timeSeriesDirNames{length(timeSeriesDirNames)+1} = tsFileName;
   end
end

% INITIALIZE MATRICES FOR STORING TIME SERIES
LHtsStore = [];
RHtsStore = [];
avgTS = [];

% NOTE STIMULUS TYPE FOR FUTURE INDEXING
stimTypeArr = [];
runOrder = '';
timeSeriesStore = [];

% LOOK AT EACH FILE IN THE TIME SERIES FOLDER
for i = 1:length(currentTimeSeriesFolder)
    % CURRENT TIME SERIES FILE
    currentTSfileName = char(currentTimeSeriesFolder(i));
    % TAKE NOTE OF CONE STIMULUS TYPE
    if strfind(currentTSfileName,'LightFlux')
        stimTypeArr(length(stimTypeArr)+1) = 1;
    elseif strfind(currentTSfileName,'L_minus_M')
        stimTypeArr(length(stimTypeArr)+1) = 2;
    elseif strfind(currentTSfileName,'_S_')
        stimTypeArr(length(stimTypeArr)+1) = 3;
    else
       stimType = []; 
    end
    
    % TAKE NOTE OF THE RUN ORDER: A OR B
    if strfind(currentTSfileName,'_A_')
        runOrder(length(runOrder)+1) = 'A';
    elseif strfind(currentTSfileName,'_B_')
        runOrder(length(runOrder)+1) = 'B';
    else
       runOrderJunk = []; 
    end
    % FIND ALL FILES CONTAINING THE FILE NAME WE WANT, AS DETERMINED BY THE
    % README FILE--GET THEIR LOCATIONS IN THE FOLDER
    tsFilesLHRH = strfind(timeSeriesDirNames,currentTSfileName);
    locationsInTSfolder = find(~cellfun(@isempty,tsFilesLHRH));
%    display(num2str(length(locationsInTSfolder)));
    % LOAD THE LEFT HEMISPHERE DATA, THEN THE RIGHT HEMISPHERE
    LHtsStruct = load([dirPathTimeSeries char(timeSeriesDirNames(locationsInTSfolder(1)))]);
    RHtsStruct = load([dirPathTimeSeries char(timeSeriesDirNames(locationsInTSfolder(2)))]);
    LHts = LHtsStruct.avgTC;
    RHts = RHtsStruct.avgTC;
    % STORE FOR PLOTTING LATER
    LHtsStore(i,:) = LHts;
    RHtsStore(i,:) = RHts;
    % MEAN OF LEFT AND RIGHT HEMISPHERE
    avgTS(i,:) = (LHts+RHts)./2;
    timeSeriesStore(size(timeSeriesStore,1)+1,:) = ((avgTS(i,:) - mean(avgTS(i,:)))./mean(avgTS(i,:))).*100;
%    timeSeriesStore(size(timeSeriesStore,1)+1,:) = avgTS(i,:);
end

%% LOAD AND PLOT STIMULUS STEP FUNCTIONS

% LOAD ALL CONTENTS OF STIMULUS DIRECTORY
stimDirContents = dir(dirPathStim);

% NUMBER OF STIMULUS FOLDERS
numberOfFolders = length(stimDirContents);

% INITIALIZE CELL CONTAINING ALL STIMULUS FOLDER NAMES
folderNameCell = {};

% LOOP OVER NUMBER OF STIMULUS FOLDERS, AND CREATE CELL WITH ALL THEIR
% NAMES
for i = 1:numberOfFolders
   miniFolderName = stimDirContents(i).name;
   if length(miniFolderName)>4 & strcmp(miniFolderName(1:4),'HERO');
       folderNameCell{length(folderNameCell)+1} = miniFolderName;
   end
end

% INITIALIZE MATRIX FOR STORING BETA VALUES
betaMatrix = [];

% STORE STIMULUS ORDER FOR A AND B
stimValuesSorted_A = [];

stimValuesSorted_B = [];

% LOOP OVER STIMULUS FOLDER NAMES
for i = 1:length(folderNameCell)
   % LOOK IN EACH RUN'S FOLDER 
   currentDirPath = [dirPathStim char(folderNameCell(i))]; 
   % GET ALL THEIR CONTENTS
   runFiles = dir(currentDirPath);
   
   % INITIALIZE MATRICES FOR STORING START TIMES AND STIMULUS VALUES
   startTimes = [];
   stimValues =[];
   
   attnFunctionRes = 10;
   
   % LOOP OVER FILES IN EACH FOLDER
   for j = 1:length(runFiles)
       lengthOfCurFile = length(runFiles(j).name);
       curFile = runFiles(j).name;
       
       % WE ARE ONLY INTERESTED IN HZ_ALL FILES
       if length(curFile)>10 & strcmp(curFile(length(curFile)-9:length(curFile)),'Hz_all.txt')
          % EXTRACT THE TEMPORAL FREQ OF THE STIMULUS FROM THE FILE NAME
          stimFile = load([currentDirPath '/' curFile]); 
          freqValueTxt = curFile(length(curFile)-11:length(curFile)-10);
          
          % SOMETIMES IT IS A TWO-DIGIT NUMBER, AND OTHER TIMES A ONE-DIGIT
          % NUMBER (0,2,4,8,16,32,OR 64)
          if str2num(freqValueTxt)
              freqValueNum = str2num(freqValueTxt);
              
          else
              freqValueNum = str2num(freqValueTxt(2));            
          end
          
          % GRAB ALL VALUES IN THE FIRST COLUMN: THESE ARE STARTING TIMES
          curTimeValue = stimFile(:,1);
          % COLLECT ALL START TIMES AND THEIR CORRESPONDING STIMULUS VALUES
          startTimes(length(startTimes)+1:length(startTimes)+length(curTimeValue)) = curTimeValue;
          stimValues(length(stimValues)+1:length(stimValues)+length(curTimeValue)) = freqValueNum;         
          
          % IF THE FILE CONTAINS ATTENTION TASK DATA
       elseif length(curFile)>20 & strcmp(curFile(length(curFile)-16:length(curFile)),'attentionTask.txt')
           attnFile = load([currentDirPath '/' curFile]);
           % RUN FIR ON THE ATTENTION START TIMES
           [hrf,~] = attentionFIR(T_R.*(1:length(avgTS(i,:))),((avgTS(i,:)-mean(avgTS(i,:)))./mean(avgTS(i,:))).*100,attnFile(:,1),lengthAttnHRF,T_R);
           hrfStore(i,:) = hrf;
           
           attnTimeValues = [];
           attnStimValues = [];
           
           % ANCHOR THE ATTENTION STEP FUNCTION AT 0
           attnStartTimes1 = attnFile(1,1);
           
           if abs(attnStartTimes1) > 0.01
              attnTimeValues(1) = 0;
              attnStimValues(1) = 0;
           end
           
           % ATTENTION TASK IS JUST A DIRAC DELTA FUNCTION TYPE DEAL
           for k = 1:size(attnFile,1)
                curTimeValue = attnFile(k,1);
                attnTimeValues = [attnTimeValues [curTimeValue linspace(curTimeValue+1e-10,curTimeValue+attnTime-1e-3,attnFunctionRes) ...
                                 curTimeValue+attnTime-1e-7]]; 
                attnStimValues = [attnStimValues [0 repmat(attnCode,[1 attnFunctionRes]) 0]];
           end
       end
       
   end
   
   % SORT THE BIG VECTOR OF START TIMES
   [startTimesSorted, stmsInd] = sort(startTimes);
   % SORT THE CORRESPONDING STIMULUS VALUES
   stimValuesSorted = stimValues(stmsInd);
   
   % INITIALIZE MATRICES FOR STORING ACTUAL STIMULUS PLOTS
   stimPlotsTimeSamples = [];
   % MARKS WHAT STIMULUS THERE WAS AT EACH TIME POINT
   stimPlotsValues = [];
   % STORES THE COSINE-WINDOWED STIMULUS
   cosWindowedStim = [];
   
   % MAKE SURE PLOT STARTS AT 0--SO CONVOLUTION WON'T RETURN NANS
   if abs(startTimesSorted(1)) > 0.001;
      stimValuesSorted = [0 stimValuesSorted];
      startTimesSorted = [0 startTimesSorted];
   end
   
   % STORES THE A AND B SEQUENCES--FOR LABELING TIME SERIES BY STIMULUS
   % PERIOD LATER
   if strfind(char(currentTimeSeriesFolder(i)),'_A_') & isempty(stimValuesSorted_A)
      startTimesSorted_A = startTimesSorted;   
      stimValuesSorted_A = stimValuesSorted;
   elseif strfind(char(currentTimeSeriesFolder(i)),'_B_') & isempty(stimValuesSorted_B)
      startTimesSorted_B = startTimesSorted;
      stimValuesSorted_B = stimValuesSorted;
   else
      stimOrderMarker = []; 
   end
   
   % HOW MANY POINTS TO SAMPLE BETWEEN STIMULUS ONSET AND OFFSET
   stepFunctionRes = 50;
   
   % DURATION OF COSINE RAMP IN SECONDS
   cosRamp = 3;
   
   % LOOP OVER ALL THE START TIMES
   for j = 1:length(startTimesSorted)
       % GRAB EACH INDIVIDUAL START TIME
       curTimeValue2 = startTimesSorted(j);
       curStimValue2 = stimValuesSorted(j);
       % MAKE 'BOX'--ADD TINY OFFSET TO MAKE SURE INTERPOLATION WORKS
       % PROPERLY
       timeValues = [curTimeValue2 linspace(curTimeValue2+1e-7,curTimeValue2+stimTime-1e-5,stepFunctionRes) ...
                     curTimeValue2+stimTime-1e-7]; 
       stimValues = [-1 repmat(curStimValue2,[1 length(timeValues)-2]) -1];
       % SHIFT TIME VALUES SO THAT STIMULUS ONSET BECOMES TIME 0
       timeForCos = timeValues - curTimeValue2;
       % REMOVE STIMULUS ONSET MARKERS FOR NOW
       timeForCos = timeForCos(2:length(timeForCos)-1);
       % DEFINE POSITIONS FOR THE LEFT RAMP AND RIGHT RAMP
       leftRampInd = timeForCos <= cosRamp;
       rightRampInd = timeForCos >= stimTime - cosRamp;
       % CREATE THE HALF COSINE
       halfCosine = ones([1 length(timeForCos)]);
       halfCosine(leftRampInd) = fliplr((cos(linspace(0,pi,sum(leftRampInd)))+1)/2);
       halfCosine(rightRampInd) = (cos(linspace(0,pi,sum(rightRampInd)))+1)/2;
       % PUT IN POINTS CORRESPONDING TO STIMULUS START AND END
       halfCosine = [-1 halfCosine -1];
       % STICK ALL THE ABOVE IN MATRICES DEFINED BEFORE THE LOOP
       stimPlotsTimeSamples = [stimPlotsTimeSamples timeValues];
       stimPlotsValues = [stimPlotsValues stimValues];
       cosWindowedStim = [cosWindowedStim halfCosine];
   end
   
   % ALL POSSIBLE STIMULI (SANS ATTENTION TASK)
   stimHz = [2 4 8 16 32 64];
   
   % INITIALIZE MATRIX OF REGRESSORS
   regMatrix = [];
   
   % MAKE HRF (TAKEN FROM GEOFF'S WINAWER MODEL CODE)
       % TOTAL DURATION IS SIMPLY LARGEST TIME VALUE
       modelDuration=floor(max(stimPlotsTimeSamples));
       modelResolution=20; 
       % TIME SAMPLES TO INTERPOLATE
       t = linspace(1,modelDuration,modelDuration.*modelResolution);
       
   if bCanonicalHRF == 1     
       % DOUBLE GAMMA HRF
       BOLDHRF = gampdf(t, param.gamma1, 1) - ...
       gampdf(t, param.gamma2, 1)/param.gammaScale;
       % scale to unit sum to preserve amplitude of y following convolution
       BOLDHRF = BOLDHRF/sum(BOLDHRF);
   else
       load([subj_name '_HRF']);
       BOLDHRF_unInterp = zeros([1 length(avgTS(i,:))]);
       hrf = hrf-hrf(1);
       BOLDHRF_unInterp(1:length(hrf)) = hrf;
       BOLDHRF = interp1(TS_timeSamples,BOLDHRF_unInterp,t);
       BOLDHRF(isnan(BOLDHRF)) = 0;
       
   end
    
   % LOOP OVER POSSIBLE STIMULUS VALUES 
   for j = 1:length(stimHz)
      % GET ALL POSITIONS WITH A GIVEN STIMULUS VALUE
      stimPositions = stimPlotsValues == stimHz(j);
      % PUT IN COSINE-WINDOWED STIMULI
      stimPreUps = zeros([1 length(stimPositions)]);
      stimPreUps(stimPositions) = cosWindowedStim(stimPositions);
      % SAMPLE AT POINTS t
      stimulusUpsampled = interp1(stimPlotsTimeSamples,stimPreUps,t,'linear','extrap');
      stimulusUpsampled = stimulusUpsampled(1:length(t));
      % CONVOLVE STIMULUS WITH HRF TO GET REGRESSOR
      regressorPreCut = conv(stimulusUpsampled,BOLDHRF);
      % CUT OFF THE EXTRA CONV VALUES--NEED TO LOOK MORE INTO THIS. CONV IS
      % WEIRD IN MATLAB
      regressor = regressorPreCut(1:length(stimulusUpsampled));
      regressorDownSampled = interp1(t,regressor,TS_timeSamples,'linear','extrap');
      % STORE THE REGRESSOR
      regMatrix(:,j) = regressorDownSampled'-mean(regressorDownSampled); 
 
%        figure;
%        plot(stimPlotsTimeSamples,stimPositions); hold on
%        plot(t,stimulusUpsampled); 
%        xlabel('time/s');
%        plot(t,regressor); title(['Stimulus and BOLD signal for ' num2str(stimHz(j)) ' Hz' ' flicker']);
%        pause;
%        close;
   end
   
   % ADD A COVARIATE OF ONES TO THE END
   regMatrix = [ones([size(regMatrix,1) 1]) regMatrix];
   
   % GET THE STEP FUNCTION FOR THE ATTENTION TASK
   attnPositions = attnStimValues == attnCode;
   attnPositions = double(attnPositions);
   % SAMPLE IT EVENLY
   attnBox = interp1(attnTimeValues,attnStimValues,t);
   % REMOVE NANS
   attnBox(isnan(attnBox)) = 0;
   % CONVOLVE WITH HRF
   attnCovariate = conv(attnBox,BOLDHRF);
   % SAMPLE TO BE THE SAME LENGTH AS OTHER REGRESSORS
   attnCovariate = attnCovariate(1:length(t));
   attnCovariateDownSampled = interp1(t,attnCovariate,TS_timeSamples,'linear','extrap');
   
   % ADD TO THE DESIGN MATRIX
   regMatrix(:,size(regMatrix,2)+1) = attnCovariateDownSampled - mean(attnCovariateDownSampled);
  
   % OBTAIN BETA WEIGHTS AND PLOT
   betaWeights = regMatrix\avgTS(i,:)'; 
   
   % BETA WEIGHTS SANS WEIGHT FOR THE FIRST REGRESSOR
   betaMatrix(i,:) = betaWeights(2:length(betaWeights)-1)./mean(avgTS(i,:));
   
   % RECONSTRUCT THE TIME SERIES ACCORDING TO MODEL, CONVERT TO %
   reconstructedTS = sum(repmat(betaWeights',[size(regMatrix,1) 1]).*regMatrix,2);
   reconstructedTS = ((reconstructedTS - mean(reconstructedTS))./mean(reconstructedTS)).*100;
   
   % STORE ALL RECONSTRUCTED TIME SERIES
   reconstructedTSmat(i,:) = reconstructedTS;

end

% SELF-EXPLANATORY VARIABLE NAMES
numberOfRuns = 12;
numRunsPerStimOrder = 6;

% CREATE CELLS FOR LABELLING PLOTS
stimValuesMatSorted_A_cell = {};
for i = 1:length(stimValuesSorted_A)
   stimValuesMatSorted_A_cell{i} = num2str(stimValuesSorted_A(i)); 
end

stimValuesMatSorted_B_cell = {};
for i = 1:length(stimValuesSorted_B)
   stimValuesMatSorted_B_cell{i} = num2str(stimValuesSorted_B(i)); 
end

%% AVERAGING A BUNCH OF THINGS TOGETHER

% CONVERT MEAN-SUBTRACTED BETA VALUES TO PERCENTAGES
LightFluxBeta =  mean(betaMatrix(stimTypeArr == 1,:)).*100;
L_minus_M_Beta = mean(betaMatrix(stimTypeArr == 2,:)).*100;
S_Beta =         mean(betaMatrix(stimTypeArr == 3,:)).*100;

% COMPUTE STANDARD ERROR
LightFluxBetaSE =  ((std(betaMatrix(stimTypeArr == 1,:)))./sqrt(numberOfRuns)).*100;
L_minus_M_BetaSE = ((std(betaMatrix(stimTypeArr == 2,:)))./sqrt(numberOfRuns)).*100;
S_BetaSE =         ((std(betaMatrix(stimTypeArr == 3,:)))./sqrt(numberOfRuns)).*100;

% AVERAGE TIME SERIES FOR EACH COMBINATION OF STIMULUS TYPE AND RUN ORDER
LightFluxAvgTS_A =  mean(timeSeriesStore(stimTypeArr == 1 & runOrder == 'A',:));
L_minus_M_AvgTS_A = mean(timeSeriesStore(stimTypeArr == 2 & runOrder == 'A',:));
S_AvgTS_A =         mean(timeSeriesStore(stimTypeArr == 3 & runOrder == 'A',:));

LightFluxAvgTS_B =  mean(timeSeriesStore(stimTypeArr == 1 & runOrder == 'B',:));
L_minus_M_AvgTS_B = mean(timeSeriesStore(stimTypeArr == 2 & runOrder == 'B',:));
S_AvgTS_B =         mean(timeSeriesStore(stimTypeArr == 3 & runOrder == 'B',:));

% STANDARD ERROR OF TIME SERIES
LightFluxStdTS_A =  (std(timeSeriesStore(stimTypeArr == 1 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder);
L_minus_M_StdTS_A = (std(timeSeriesStore(stimTypeArr == 2 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder);
S_StdTS_A =         (std(timeSeriesStore(stimTypeArr == 3 & runOrder == 'A',:)))./sqrt(numRunsPerStimOrder);

LightFluxStdTS_B =  (std(timeSeriesStore(stimTypeArr == 1 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder);
L_minus_M_StdTS_B = (std(timeSeriesStore(stimTypeArr == 2 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder);
S_StdTS_B =         (std(timeSeriesStore(stimTypeArr == 3 & runOrder == 'B',:)))./sqrt(numRunsPerStimOrder);

% DO THE SAME FOR THE 'RECONSTRUCTED' TIME SERIES'
LightFluxAvgTS_Model_A =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'A',:));
L_minus_M_AvgTS_Model_A = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'A',:));
S_AvgTS_Model_A =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'A',:));

LightFluxAvgTS_Model_B =  mean(reconstructedTSmat(stimTypeArr == 1 & runOrder == 'B',:));
L_minus_M_AvgTS_Model_B = mean(reconstructedTSmat(stimTypeArr == 2 & runOrder == 'B',:));
S_AvgTS_Model_B =         mean(reconstructedTSmat(stimTypeArr == 3 & runOrder == 'B',:));

yLimits = [min([LightFluxBeta L_minus_M_Beta S_Beta]) max([LightFluxBeta L_minus_M_Beta S_Beta])];

%% TTF AND HRF PLOTS
[wftd1, fp1] = fitWatsonToTTF_errorGuided(stimHz,LightFluxBeta,LightFluxBetaSE,1); hold on
errorbar(stimHz,LightFluxBeta,LightFluxBetaSE,'ko');
set(gca,'FontSize',15);
set(gca,'Xtick',stimHz);
title('Light flux');
[wftd2, fp2] = fitWatsonToTTF_errorGuided(stimHz,L_minus_M_Beta,L_minus_M_BetaSE,1); hold on
errorbar(stimHz,L_minus_M_Beta,L_minus_M_BetaSE,'ko');
set(gca,'FontSize',15);
set(gca,'Xtick',stimHz);
title('L - M');
[wftd3, fp3] = fitWatsonToTTF_errorGuided(stimHz,S_Beta,S_BetaSE,1); hold on
errorbar(stimHz,S_Beta,S_BetaSE,'ko');
set(gca,'FontSize',15);
set(gca,'Xtick',stimHz);
title('S');

figure;
errorbar(T_R.*(1:lengthAttnHRF)-1,mean(hrfStore),std(hrfStore)./sqrt(size(hrfStore,1)),'LineWidth',2);
xlabel('Time/s'); ylabel('% signal change'); title('HRF obtained using FIR');
set(gca,'FontSize',15);

%% TIME SERIES PLOTS
figure;
set(gcf,'Position',[156 372 1522 641])

subplot(3,2,1)
plotLinModelFits(T_R.*(1:length(LightFluxAvgTS_A)),LightFluxAvgTS_A,LightFluxAvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,LightFluxStdTS_A);
title('Light flux A'); xlabel('Time / s'); ylabel('% signal change');

subplot(3,2,3)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_A)),L_minus_M_AvgTS_A,L_minus_M_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,L_minus_M_StdTS_A);
title('L - M A'); xlabel('Time / s'); ylabel('% signal change');
              
subplot(3,2,5)
plotLinModelFits(T_R.*(1:length(S_AvgTS_A)),S_AvgTS_A,S_AvgTS_Model_A, ...
                 startTimesSorted_A,stimValuesMatSorted_A_cell,stimValuesSorted_A,S_StdTS_A);
title('S A'); xlabel('Time / s'); ylabel('% signal change');
 
subplot(3,2,2)
plotLinModelFits(T_R.*(1:length(LightFluxAvgTS_B)),LightFluxAvgTS_B,LightFluxAvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,LightFluxStdTS_B);
title('Light flux B');
 
subplot(3,2,4)
plotLinModelFits(T_R.*(1:length(L_minus_M_AvgTS_B)),L_minus_M_AvgTS_B,L_minus_M_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,L_minus_M_StdTS_B);
title('L - M B');
              
subplot(3,2,6)
plotLinModelFits(T_R.*(1:length(S_AvgTS_B)),S_AvgTS_B,S_AvgTS_Model_B, ...
                 startTimesSorted_B,stimValuesMatSorted_B_cell,stimValuesSorted_B,S_StdTS_B);
title('S B');