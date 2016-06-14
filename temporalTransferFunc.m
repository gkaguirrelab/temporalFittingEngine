%%

% CREATES TEMPORAL TRANSFER FUNCTIONS. CODE STILL HAS SOME NUMERIC ISSUES.
% ATTENTION TASK NOT YET ACCOUNTED FOR. 

% SPECIFY THE SUBJECT AND DATE BELOW. 

% MAKES PLOTS ONE BY ONE--PRESS ANY KEY TO CYCLE THROUGH THEM

%% Identify the user
 if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
 end


%% SPECIFY SUBJECT AND SESSION, AND DROPBOX FOLDER

subj_name = 'HERO_asb1';
%     'HERO_asb1' 
%     'HERO_gka1'


session = '041416';
%     '041416' ...
%     '041516' ...

% PATH TO LOCAL DROPBOX
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'];

%% HRF PARAMETERS (GRABBED FROM WINAWER MODEL CODE)

param = struct;
% parameters of the double-gamma hemodynamic filter (HRF)
param.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in seconds)
param.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in seconds)
param.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
        
%%

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

% SUBJECT AND DATE DETERMINE WHICH TIME SERIES FILES WE LOAD, AND THE ORDER
% THEY ARE PLOTTED IN
if strcmp(subj_name,'HERO_asb1') & strcmp(session,'041416')
    currentTimeSeriesFolder = tsFileNamesASB1_DAY1;
elseif strcmp(subj_name,'HERO_asb1') & strcmp(session,'041516')
    currentTimeSeriesFolder = tsFileNamesASB1_DAY2;
elseif strcmp(subj_name,'HERO_gka1') & strcmp(session,'041416')
    currentTimeSeriesFolder = tsFileNamesGKA1_DAY1;
elseif strcmp(subj_name,'HERO_gka1') & strcmp(session,'041516')
    currentTimeSeriesFolder = tsFileNamesGKA1_DAY2;
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
LHtsMat = [];
RHtsMat = [];

% LOOK AT EACH FILE IN THE TIME SERIES FOLDER
for i = 1:length(currentTimeSeriesFolder)
    % CURRENT TIME SERIES FILE
    currentTSfileName = char(currentTimeSeriesFolder(i));
    % FIND ALL FILES CONTAINING THE FILE NAME WE WANT, AS DETERMINED BY THE
    % README FILE--GET THEIR LOCATIONS IN THE FOLDER
    tsFilesLHRH = strfind(timeSeriesDirNames,currentTSfileName);
    locationsInTSfolder = find(~cellfun(@isempty,tsFilesLHRH));
    % LOAD THE LEFT HEMISPHERE DATA, THEN THE RIGHT HEMISPHERE
    LHtsStruct = load([dirPathTimeSeries char(timeSeriesDirNames(locationsInTSfolder(1)))]);
    RHtsStruct = load([dirPathTimeSeries char(timeSeriesDirNames(locationsInTSfolder(2)))]);
    LHts = LHtsStruct.avgTC;
    RHts = RHtsStruct.avgTC;
    % STORE FOR PLOTTING LATER
    LHtsMat(i,:) = LHts;
    RHtsMat(i,:) = RHts;
end

% LOAD ALL CONTENTS OF STIMULUS DIRECTORY
files1 = dir(dirPathStim);

% NUMBER OF STIMULUS FOLDERS
numberOfFolders = length(files1);

% INITIALIZE CELL CONTAINING ALL STIMULUS FOLDER NAMES
folderNameCell = {};

% DURATION OF STIMULUS (ALWAYS THE SAME)
stimTime = 12;

% LOOP OVER NUMBER OF STIMULUS FOLDERS, AND CREATE CELL WITH ALL THEIR
% NAMES
for i = 1:numberOfFolders
   miniFolderName = files1(i).name;
   if length(miniFolderName)>4 & strcmp(miniFolderName(1:4),'HERO');
       folderNameCell{length(folderNameCell)+1} = miniFolderName;
   end
end

% LOOP OVER STIMULUS FOLDER NAMES
for i = 1:length(folderNameCell)
   % LOOK IN EACH RUN'S FOLDER 
   currentDirPath = [dirPathStim char(folderNameCell(i))]; 
   % GET ALL THEIR CONTENTS
   runFiles = dir(currentDirPath);
   
   % INITIALIZE MATRICES FOR STORING START TIMES AND STIMULUS VALUES
   startTimesMat = [];
   stimValuesMat =[];
   
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
          startTimesMat(length(startTimesMat)+1:length(startTimesMat)+length(curTimeValue)) = curTimeValue;
          stimValuesMat(length(stimValuesMat)+1:length(stimValuesMat)+length(curTimeValue)) = freqValueNum;         
          
          % IF THE FILE CONTAINS ATTENTION TASK DATA * NEED TO DEAL WITH THIS DATA PROPERLY *
       elseif length(curFile)>20 & strcmp(curFile(length(curFile)-16:length(curFile)),'attentionTask.txt')
           attnFile = load([currentDirPath '/' curFile]); 
           
           attnTimeValues = [];
           attnStimValues = [];
           
           % ATTENTION TASK IS JUST A DIRAC DELTA FUNCTION TYPE DEAL
           for k = 1:size(attnFile,1)
                curTimeValue = attnFile(k,1);
                attnTimeValues = [attnTimeValues [curTimeValue curTimeValue+1e-10 curTimeValue+1e-8]]; 
                attnStimValues = [attnStimValues [0 96 0]];
           end
       end
       
   end
   
   % SORT THE BIG VECTOR OF START TIMES
   [startTimesMatSorted, stmsInd] = sort(startTimesMat);
   % SORT THE CORRESPONDING STIMULUS VALUES
   stimValuesMatSorted = stimValuesMat(stmsInd);
   
   % INITIALIZE MATRICES FOR STORING ACTUAL STIMULUS PLOTS
   timeValuesMatFinal = [];
   stimValuesMatFinal = [];
   
   % LOOP OVER ALL THE START TIMES
   for j = 1:length(startTimesMatSorted)
       % GRAB EACH INDIVIDUAL START TIME
       curTimeValueFinal = startTimesMatSorted(j);
       curStimValueFinal = stimValuesMatSorted(j);
       % MAKE 'BOX'--ADD TINY OFFSET TO MAKE SURE INTERPOLATION WORKS
       % PROPERLY
       timeValues = [curTimeValueFinal curTimeValueFinal+1e-7 ...
                    curTimeValueFinal+stimTime-1e-3 curTimeValueFinal-1e-7]; 
       stimValues = [0 curStimValueFinal curStimValueFinal 0];
       % STICK IN MATRICES DEFINED BEFORE THE LOOP
       timeValuesMatFinal = [timeValuesMatFinal timeValues];
       stimValuesMatFinal = [stimValuesMatFinal stimValues];
   end
   
   % MAKE SURE PLOT STARTS AT 0--SO CONVOLUTION WON'T RETURN NANS
   if timeValuesMatFinal(1) ~= 0
      stimValuesMatFinal = [0 stimValuesMatFinal];
      timeValuesMatFinal = [0 timeValuesMatFinal];
   end
   
   % ALL POSSIBLE STIMULI (SANS ATTENTION TASK)
   stimHz = [0 2 4 8 16 32 64];
   
   % INITIALIZE MATRIX OF COVARIATES
   regMatrix = [];
   
   % MAKE HRF (TAKEN FROM GEOFF'S WINAWER MODEL CODE)
   % TOTAL DURATION IS SIMPLY LARGEST TIME VALUE
   modelDuration=max(timeValuesMatFinal);
   % I HAVE INFERRED THAT THE TEMPORAL RESOLUTION IS 1 FROM THE TIME SERIES, 
   % BUT IS IT ACTUALLY 0.8?
   modelResolution=1; 
   % TIME SAMPLES TO INTERPOLATE
   t = linspace(1,modelDuration,modelDuration.*modelResolution);
   % DOUBLE GAMMA HRF
    BOLDHRF = gampdf(t, param.gamma1, 1) - ...
    gampdf(t, param.gamma2, 1)/param.gammaScale;
    % scale to unit sum to preserve amplitude of y following convolution
    BOLDHRF = BOLDHRF/sum(BOLDHRF);
    
   % LOOP OVER POSSIBLE STIMULUS VALUES 
   for j = 1:length(stimHz)
      % GET ALL INSIXWA WITH A GIVEN STIMULUS VALUE
      stimPositions = stimValuesMatFinal == stimHz(j); 
      stimPositions = double(stimPositions);
      % SAMPLE AT POINTS t
      stimulusCovariate = interp1(timeValuesMatFinal,stimPositions,t);
      % CONVOLVE STIMULUS WITH HRF TO GET REGRESSOR
      regressor = conv(stimulusCovariate,BOLDHRF);
      % CUT OFF THE EXTRA CONV VALUES--NEED TO LOOK MORE INTO THIS. CONV IS
      % WEIRD IN MATLAB
      regressor = regressor(1:length(stimulusCovariate));
      % STORE THE REGRESSOR
      regMatrix(:,j) = regressor';
      % STORE THE COVARIATES
      SCmat(:,j) = stimulusCovariate;
   end
   
   % OBTAIN BETA WEIGHTS AND PLOT
   betaWeights = regMatrix\LHtsMat(i,:)';
   
   figure;
   plot(stimHz,betaWeights);
   pause;
   close;
%    figure;
%    set(gcf,'Position',[439 222 1029 876]);
%    subplot(3,1,3);
%    plot(timeValuesMatFinal,stimValuesMatFinal); hold on
%    plot(attnTimeValues,attnStimValues);
%    xlabel('Time(s)'); ylabel('Stimulus frequency (Hz)');
%    title('Stimulus');
%    set(gca,'FontSize',15);
%    subplot(3,1,2);
%    plot(LHtsMat(i,:));
%    title('Left hemisphere BOLD response');
%    set(gca,'FontSize',15);
%    subplot(3,1,1);
%    plot(RHtsMat(i,:));
%    title('Right hemisphere BOLD response');
%    set(gca,'FontSize',15);
%    pause;
%    close;
end