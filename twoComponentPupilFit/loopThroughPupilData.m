


% Identify the user
if isunix
    [~, user_name] = system('whoami') ;
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%') ;
end

% Path to local Dropbox
localDropboxDir = ['/Users/',strtrim(user_name),'/Dropbox-Aguirre-Brainard-Lab/'] ;

% Create the output directory
outPathStem=[localDropboxDir 'MELA_analysis/MelanopsinMR/pupilModelFits/'];
mkdir(outPathStem);

% subjects and stimuli
subjectNames={'HERO_asb1';'HERO_gka1';'HERO_aso1';'HERO_mxs1'};
modulationDirections={'MELA_analysis/MelanopsinMRMaxLMSCRF';'MELA_analysis/MelanopsinMRMaxMelCRF'};
contrastLabels={'25%';'50%';'100%';'200%';'400%'};
modulationLabels={'LMS';'Mel'};

% labels for the parameters of the model
paramLabels={'responseOnset','gammaShape','sustainAmplitude',...
    'sustainTau','persistentAmplitude','persistentT50',...
    'persistentAlpha'};


%% Obtain the pupil stimulus vector
stimulus = createPupilStimulusVector();

% Index to count up the results
resultsIndex=1;

% Loop through stimuli
for stim=1:length(modulationLabels)
    
    % Loop through the subjects
    for sub=1:length(subjectNames)
        
        % Define a path to Time Series data
        filePathTimeSeriesDir = [localDropboxDir char(modulationDirections(stim)) '/' ...
            char(subjectNames(sub)) '/'] ;
        filePathTimeSeries=dir([filePathTimeSeriesDir '*MeanTimeSeries.csv']);
        filePathTimeSeries=[filePathTimeSeriesDir filePathTimeSeries(1).name];
        
        % Load the data file and extract the time vector
        rawDataFile=load(filePathTimeSeries);
        t=rawDataFile(:,1); % first column has time point vector
        
        for contrast=1:length(contrastLabels)
            data=rawDataFile(:,contrast+1); % need to skip the first column, which has time data
            stdErrorVector=rawDataFile(:,contrast+1+5); % need to skip the first column, which has time data
            
            [yPupil, fitParams, fitFigureHandle]=fitPupilModelToData(t,data,stdErrorVector,stimulus);
            
            % Save the figure and fitParams
            figureFileNameStem=[char(subjectNames(sub)) '_' char(modulationLabels(stim)) '_' char(contrastLabels(contrast))];
            saveas(fitFigureHandle,[outPathStem figureFileNameStem],'png');
            fitParamResults(resultsIndex,:)=fitParams;
            
            % Calculate the R2 of the model
            fitVarianceExplained(resultsIndex,1)=corr2(data,yPupil)^2;
            
            % Make some dumb vectors to simplify later data presentation
            stimVec(resultsIndex,1)=stim;
            subVec(resultsIndex,1)=sub;
            conVec(resultsIndex,1)=contrast;
            
            % Increment the results index
            resultsIndex=resultsIndex+1;
            
            % close the figure
            close(fitFigureHandle);
            
            % assemble an average pupil response across contrast levels
            if contrast==1
                averageData=data;
                averageStdErrorVector=stdErrorVector;
            else
                averageData=averageData+data;
                averageStdErrorVector=averageStdErrorVector+stdErrorVector;
            end % checking if we are on first contrast loop
        end % contrast levels
        
        clear data
        clear stdErrorVector
        
        % Plot and fit the average pupil response
        averageData=averageData/length(contrastLabels);
        averageStdErrorVector=averageStdErrorVector/length(contrastLabels);
        [yPupil, fitParams, fitFigureHandle]=fitPupilModelToData(t,averageData,averageStdErrorVector,stimulus);
        
        % Save the figure and fitParams
        figureFileNameStem=[char(subjectNames(sub)) '_' char(modulationLabels(stim)) '_AcrossContrastAverage' ];
        saveas(fitFigureHandle,[outPathStem figureFileNameStem],'png');
        fitParamAverageVector(sub,:)=fitParams;
        
    end % subjects
    
end % stimuli


% Save the variables that contain the fit params and variance explained
variableFileName=['fitParamResults.mat'];
save([outPathStem variableFileName],'fitParamResults','fitVarianceExplained','stimVec','subVec','conVec','fitParamAverageVector');

%% Make some summary plots

% Colors for the plots
modulationPlotColors={'k','r'};

% Clear out residual holds from the prior plots
close all
hold off

% Plot each of the parameters vs. contrast
x=linspace(1,length(contrastLabels),length(contrastLabels));
for param=1:7
    paramFigureHandle=figure;
    hold on
    clear r;
    for stim=1:length(modulationLabels)
        clear meanParam;
        for contrast=1:length(contrastLabels)
            index=find((conVec==contrast).*(stimVec==stim));
            meanParam(contrast)=mean(fitParamResults(index,param));
            semParam(contrast)=std(fitParamResults(index,param))/sqrt(length(index));
        end % loop across contrast levels
        r(stim)=errorbar(x,meanParam,semParam,char(modulationPlotColors(stim)),'LineWidth',1,'LineStyle','none');
        plot(x,meanParam,'Color',char(modulationPlotColors(stim)),'Marker','o','MarkerSize',5);
        xlim([0,length(contrastLabels)+1]);
        ax = gca;
        ax.XTickLabels=[' ';contrastLabels;' '];
    end % stim
    title(paramLabels(param));
    xlabel('Contrast');
    ylabel('Mean ±SEM param value');
    legend(r, modulationLabels); legend boxoff;
    % save the figure
    figureFileNameStem=['meanParam_' char(paramLabels(param))];
    saveas(paramFigureHandle,[outPathStem figureFileNameStem],'png');
    hold off
end % param

% Plot the ratio of sustained vs. persistent amplitudes across contrast,
% and do so for Mel and LMS.

paramFigureHandle=figure;
for stim=1:length(modulationLabels)
    clear xAxisParam;
    clear xAxisSEM;
    paramIndex=3; % This is the amplitude of the sustained component
    for contrast=1:length(contrastLabels)
        index=find((conVec==contrast).*(stimVec==stim));
        xAxisParam(contrast)=mean(fitParamResults(index,paramIndex));
        xAxisSEM(contrast)=std(fitParamResults(index,paramIndex))/sqrt(length(index));
    end % loop across contrast levels
    clear yAxisParam;
    clear yAxisSEM;
    paramIndex=5; % This is the amplitude of the persistent component
    for contrast=1:length(contrastLabels)
        index=find((conVec==contrast).*(stimVec==stim));
        yAxisParam(contrast)=mean(fitParamResults(index,paramIndex));
        yAxisSEM(contrast)=std(fitParamResults(index,paramIndex))/sqrt(length(index));
    end % loop across contrast levels
    % plot one parameter vs. the other
    r(stim)=plot(xAxisParam,yAxisParam,'Color','none','MarkerFaceColor',char(modulationPlotColors(stim)),'Marker','o','MarkerSize',10);
    hold on
    herrorbar(xAxisParam,yAxisParam,xAxisSEM,['.' char(modulationPlotColors(stim))]);
    errorbar(xAxisParam,yAxisParam,yAxisSEM,['.' char(modulationPlotColors(stim))]);
    % add a 2 dimensional polynomial fit
    fitLineYAxis=polyval(polyfit(xAxisParam,yAxisParam,2),linspace(min(xAxisParam),max(xAxisParam),100));
    plot(linspace(min(xAxisParam),max(xAxisParam),100),fitLineYAxis,'-b');
end % stim
title('Ratio sustained amp to persistent amp [±SEM]');
xlabel('sustained amp');
ylabel('persistent amp');
legend(r, modulationLabels);
hold off;
% save the figure
figureFileNameStem=['paramPlot_sustainedAmp_v_persistentAmp'];
saveas(paramFigureHandle,[outPathStem figureFileNameStem],'png');

% Plot the ratio of delay at onset vs. gamma shape across contrast,
% and do so for Mel and LMS. This captures how "sluggish" the initial
% response is.
paramFigureHandle=figure;
for stim=1:length(modulationLabels)
    clear xAxisParam;
    clear xAxisSEM;
    paramIndex=1; % Latency of response onset
    for contrast=1:length(contrastLabels)
        index=find((conVec==contrast).*(stimVec==stim));
        xAxisParam(contrast)=mean(fitParamResults(index,paramIndex));
        xAxisSEM(contrast)=std(fitParamResults(index,paramIndex))/sqrt(length(index));
    end % loop across contrast levels
    clear yAxisParam;
    clear yAxisSEM;
    paramIndex=2; % Shape parameter of the gamma function
    for contrast=1:length(contrastLabels)
        index=find((conVec==contrast).*(stimVec==stim));
        yAxisParam(contrast)=mean(fitParamResults(index,paramIndex));
        yAxisSEM(contrast)=std(fitParamResults(index,paramIndex))/sqrt(length(index));
    end % loop across contrast levels
    r(stim)=plot(xAxisParam,yAxisParam,'Color','none','MarkerFaceColor',char(modulationPlotColors(stim)),'Marker','o','MarkerSize',10);
    hold on
    herrorbar(xAxisParam,yAxisParam,xAxisSEM,['.' char(modulationPlotColors(stim))]);
    errorbar(xAxisParam,yAxisParam,yAxisSEM,['.' char(modulationPlotColors(stim))]);
end % stim
title('Ratio latency to gamma shape [±SEM]');
xlabel('latency of onset');
ylabel('gamma shape parameter');
legend(r, modulationLabels);
hold off;
% save the figure
figureFileNameStem=['paramPlot_latency_v_gammaShape'];
saveas(paramFigureHandle,[outPathStem figureFileNameStem],'png');
