% This Script creates Temporal Transer Functions 
% using a linear model and FIR with cosine windows.

%% Variable name legend 

% ts        = time series
% LH & RH   = left & right hemisphere
% stim      = stimulus
% AVG       = average
% Arr       = array (superfluous--one of Ben's coding habit)
% Mat       = matrix
% cur       = current
% position  = index

%% Identify the user
 if isunix
    [~, user_name] = system('whoami') ; % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%') ; % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
 end


%% Specify Subject & Session, With Dropbox Folder

subj_name = 'HERO_gka1' ; 
% *** Subject Pool ***
%     'HERO_asb1' 
%     'HERO_gka1'


session = 'all' ;
% *** Dates ***
%     '041416' ...
%     '041516' ...

% Path to local Dropbox
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'] ;

% Define path to folder for one sunject on one date
dirPathStim = [localDropboxDir 'MELA_analysis/HCLV_Photo_7T/mriTemporalFitting_data/' ...
           subj_name '/' session '/' 'Stimuli/'] ;
        
%% Defining Paths, Order, etc

% load time series
[avgTS, avgTSprc, tsFileNames, stimTypeArr, runOrder] = loadTimeSeriesData(subj_name,session);

% get all stimulus values and start times, as well as the attention task
% start times
[startTimesSorted, stimValuesSorted, attnStartTimes] = orderStartTimes(subj_name,session);
