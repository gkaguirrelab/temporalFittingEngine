function [ fitParamsMean, fitParamsSD ] = t_IAMPDesignFMRIExperiment(varargin)
% function [ paramsFit ] = t_IAMPDesignFMRIExperiment(varargin)
%
% This function explores the how different orderings of events in time
%  in an fMRI experiment yield better or worse estimates of response
%  parameters.
%
% Optional key/value pairs
%	generatePlots - true/fale (default true).  Make plots?
%	deltaT - temporal resolution of the simulation (in msecs)
%   totalTime - total duration of the simulation (in msecs)
%   trialTime - duration of each stimulus presentation (in msecs)
%   stimResponseVec - a column vector of amplitudes of neural response to
%       be modeled for each stimulus type
%   stimLabels - cell array of labels for each stimulus
%   hrfParams - a vector of params that define the double gamma HRF model,
%       corresponding to gamma1, gamma2, and gammaScale
%   noiseSd - the standard deviation of the noise
%   noiseInverseFrequencyPower - the exponent alpha that defines the
%       autocorrelated noise present in the fMRI data, following the
%       the equation: noise = 1/|freq|^alpha
%   nSimulations - the number of simulations to conduct to obtain
%       confidence intervals on the fit params

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.addParameter('deltaT',100,@isnumeric);
p.addParameter('totalTime',330000,@isnumeric);
p.addParameter('trialTime',12000,@isnumeric);
p.addParameter('stimResponseVec',[0.001; .25; .5; .75; 1],@isnumeric);
p.addParameter('stimLabels',{'0%' '25%' '50%' '100%' '200%'},@iscell);
p.addParameter('hrfParams',[6,12,10],@isnumeric);
p.addParameter('noiseSd',0.4,@isnumeric);
p.addParameter('noiseInverseFrequencyPower',1.0,@isnumeric);
p.addParameter('nSimulations',100,@isnumeric);
p.parse(varargin{:});

%% Derive values from and sanity check the input

if length(p.Results.stimResponseVec) ~= length(p.Results.stimLabels)
    error('The length of the response vec and the stimulus labels is not the same');
end
nStimTypes = length(p.Results.stimResponseVec);

nTrials = floor(p.Results.totalTime / p.Results.trialTime);
if nTrials < nStimTypes
    error('The design does not allow enough time to present all stimulus types');
end

%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

%% Temporal domain of the stimulus
stimulusStruct.timebase = linspace(0,p.Results.totalTime-p.Results.deltaT,p.Results.totalTime/p.Results.deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = p.Results.hrfParams(1);   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = p.Results.hrfParams(2);  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = p.Results.hrfParams(3); % scaling factor between the positive and negative gamma componenets

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
kernelStruct.timebase=linspace(0,15999,16000);
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% Normalize the kernel to have unit amplitude
[ kernelStruct ] = normalizeKernelArea( kernelStruct );

%% Set up the simulated response parameters
defaultParamsInfo.nInstances = nStimTypes;
params0 = temporalFit.defaultParams('defaultParamsInfo', defaultParamsInfo);

% Populate the amplitude params with the passed values
params0.paramMainMatrix=p.Results.stimResponseVec;
params0.noiseSd = p.Results.noiseSd;
params0.noiseInverseFrequencyPower = p.Results.noiseInverseFrequencyPower;


%% Loop through simulations
for ss = 1:p.Results.nSimulations
    
    % Create a random ordering of events of the specified duration
    notDoneFlag = true;
    % Keep drawing stimulus orders until we have at least one of every type
    while notDoneFlag
        stimTypes = datasample(1:1:nStimTypes,nTrials,'replace',true)';
        if length(unique(stimTypes)) == nStimTypes
            notDoneFlag = false;
        end
    end
    
    stimulusStruct.meta.stimLabels = p.Results.stimLabels;
    stimulusStruct.meta.stimTypes = stimTypes;
    
    for ii=1:nStimTypes
        stimulusStruct.values(ii,:)=zeros(1,nTimeSamples);
        eventTimes = (find(stimTypes == ii)-1)*p.Results.trialTime;
        for ee=1:length(eventTimes)
            stimulusStruct.values(ii,floor(eventTimes(ee)/p.Results.deltaT)+1:floor(eventTimes(ee)/p.Results.deltaT+floor(p.Results.trialTime/p.Results.deltaT)))=1;
        end
    end

    % Create the simulated response
    simulatedResponseStruct = temporalFit.computeResponse(params0,stimulusStruct,kernelStruct,'AddNoise',true);
    
    % Construct a packet and model params
    thePacket.stimulus = stimulusStruct;
    thePacket.response = simulatedResponseStruct;
    thePacket.kernel = kernelStruct;
    thePacket.metaData = [];
    
    % Obtain the fit params
    [paramsFit,~,modelResponseStruct] = ...
        temporalFit.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo, ...
        'searchMethod','linearRegression');
    
    paramValsAcrossSimulations(ss,:) = paramsFit.paramMainMatrix';
    
end % loop over simulations

% Obtain the mean and SD of the parameters across simulations
fitParamsMean = mean(paramValsAcrossSimulations);
fitParamsSD = std(paramValsAcrossSimulations);

% Plot the last of the fit results
if p.Results.generatePlots
    temporalFit.plot(simulatedResponseStruct,'Color',[0 1 0],'NewWindow',true,'DisplayName','last simulated signal');
    temporalFit.plot(modelResponseStruct,'Color',[1 0 0],'NewWindow',false,'DisplayName','last model fit');
    legend('show');legend('boxoff');
    hold off
end

% Plot of simulated vs. recovered parameter values
if p.Results.generatePlots
    figure
    errorbar(params0.paramMainMatrix,fitParamsMean,fitParamsSD,'.r','markerfacecolor',[1 0 0])
    xlabel('simulated stimulus neural amplitudes') % x-axis label
    ylabel('estimated stimulus neural amplitudes') % y-axis label
    ylim([ (min(fitParamsMean) - max(fitParamsSD)) (max(fitParamsMean) + max(fitParamsSD)) ]);
    xlim([min(p.Results.stimResponseVec)-.5*min(diff(p.Results.stimResponseVec)) max(p.Results.stimResponseVec)+.5*min(diff(p.Results.stimResponseVec))]);
end

end % function