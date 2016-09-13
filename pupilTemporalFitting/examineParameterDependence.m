% Examine parameter dependence on response scaling
%
% Ideally, the relative amplitude of the sustained and persistent
% components should be unchanged with scaling of the overal pupil response,
% and the time-constant values should be unaltered.
%
% This routine conducts fitting on a pupil response that has been scaled
% between 20% and 200% of the original size. We then plot the ratio of the
% amplitudes of the two response components, as well as the time constants,
% across scaling choices. The resulting plots should have zero slope.
%

% Identify the user
if isunix
    [~, user_name] = system('whoami') ;
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%') ;
end

% Path to local Dropbox
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'] ;

% subjects and stimuli
subjectNames={'HERO_asb1';'HERO_gka1';'HERO_aso1';'HERO_mxs1'};
modulationDirections={'MELA_analysis/MelanopsinMRMaxLMSCRF';'MELA_analysis/MelanopsinMRMaxMelCRF'};
contrastLabels={'25%';'50%';'100%';'200%';'400%'};
modulationLabels={'LMS';'Mel'};

%% Obtain the pupil stimulus vector
stimulus = createPupilStimulusVector();

loopIndex=1;

% Loop through stimuli
for stim=1:1   %length(modulationLabels)
    
    % Loop through the subjects
    for sub=4:4  %length(subjectNames)
        
        % Define a path to Time Series data
        filePathTimeSeriesDir = [localDropboxDir char(modulationDirections(stim)) '/' ...
            char(subjectNames(sub)) '/'] ;
        filePathTimeSeries=dir([filePathTimeSeriesDir '*MeanTimeSeries.csv']);
        filePathTimeSeries=[filePathTimeSeriesDir filePathTimeSeries(1).name];
        
        % Load the data file and extract the time vector
        rawDataFile=load(filePathTimeSeries);
        t=rawDataFile(:,1); % first column has time point vector
        
        for contrast=5:5%length(contrastLabels)
            data=rawDataFile(:,contrast+1); % need to skip the first column, which has time data
            stdErrorVector=rawDataFile(:,contrast+1+5);
            
            for scales=1:10
                scaleFactor=scales/5;
                [~, fitParams, ~] = fitPupilModelToData(t, data*scaleFactor, stdErrorVector*scaleFactor, stimulus);
                amplitudeExponential(loopIndex,scales)=fitParams(3);
                amplitudeSaturation(loopIndex,scales)=fitParams(5);
            end % scale factors
            loopIndex=loopIndex+1;
        end % contrast levels
        
    end % subjects
    
end % stimuli

hold off
figure
plot(amplitudeExponential./amplitudeSaturation);
