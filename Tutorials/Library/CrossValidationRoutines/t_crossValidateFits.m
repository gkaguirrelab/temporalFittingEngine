function [ xValFitStructure ] = t_crossValidateFits(varargin)
% function [ xValFitStructure ] = t_crossValidateFits()
%
% Demonstrate cross-validation of simulated data using the BTRM model
%
% Optional key/value pairs
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});

%% Set the random number generator to default
rng('default');

%% Construct the model object
tfeHandle = tfeIAMP('verbosity','none');

%% Temporal domain of the stimulus
deltaT = 100; % in msecs
totalTime = 120000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus values.
% We create here an impulse of stimulation. There will be ten instances.
nInstances=9;
eventDuration=100;
ISIduration=12000; % msecs

defaultParamsInfo.nInstances=nInstances;

for ii=1:nInstances
    stimulusStruct.values(ii,:)=zeros(1,nTimeSamples);
    stimulusStruct.values(ii,(ii*ISIduration)/deltaT+1:(ii*ISIduration)/deltaT+eventDuration)=1;
    
end % loop over instances building the stimulus values

%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 4;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
kernelStruct.timebase=stimulusStruct.timebase;

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct.values=kernelStruct.values-kernelStruct.values(1);
kernelStruct=normalizeKernelArea(kernelStruct);

%% Get the default forward model parameters
params0 = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
params0.noiseSd = 0.1;

%% Create the packetCellArray

% We will now create a set of nPackets. Each packet will have two
% stimulus types, each with a different expected amplitude and tau.
stimLabels={'stimA','stimB'};
stimTypes=[1 2 1 2 1 2 1 2 1 2]';
stimMeansAmplitude=[1,0.5];
nPackets=6;

for pp=1:nPackets
    
    % ensure that the packet from the pior loop is empty
    thePacket=[];
    
    % start the packet assembly
    thePacket.stimulus = stimulusStruct;
    thePacket.kernel = kernelStruct;
    thePacket.metaData = [];
    
    % Randomize the order of the stimTypes
    ix = randperm(nInstances);
    thePacket.stimulus.metaData.stimTypes=stimTypes(ix);
    thePacket.stimulus.metaData.stimLabels=stimLabels;
    
    % Create some params to define the simulated data for this packet
    paramsLocal=params0;
    for ii=1:nInstances
        paramsLocal.paramMainMatrix(ii,1)=stimMeansAmplitude(thePacket.stimulus.metaData.stimTypes(ii));
    end
    
    % Generate the simulated response
    simulatedResponseStruct = tfeHandle.computeResponse(paramsLocal,thePacket.stimulus,thePacket.kernel,'AddNoise',true);
    
    % Add the simulated response to this packet
    thePacket.response=simulatedResponseStruct;
    
    % Add this packet to the growing cell array
    packetCellArray{pp}=thePacket;
end % loop over number of packets to be created

%% Conduct the cross validation
[ xValFitStructure, ~, ~ ] = crossValidateFits( packetCellArray, tfeHandle, ...
    'partitionMethod', 'twentyPercent', ...
    'maxPartitions' , 2, ...
    'aggregateMethod', 'mean',...
    'verbosity', 'full',...
    'errorType','rmse');

%% Report the results

fprintf(['Modeled noise: ' strtrim(num2str(params0.noiseSd)) '\n']);
fprintf(['Observed error (train): ' strtrim(num2str(mean(xValFitStructure.trainfVals))) '\n']);
fprintf(['Observed error (test): ' strtrim(num2str(mean(xValFitStructure.testfVals))) '\n']);


end % function