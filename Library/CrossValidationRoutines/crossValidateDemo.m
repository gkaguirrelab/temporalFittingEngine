% crossValidateDemo
%
% Demonstrate cross-validation of simulated data using the BTRM model
%

%% Clear and close
clear; close all;

%% Construct the model object
tfeHandle = tfeBTRM('verbosity','none');

%% Temporal domain of the stimulus
deltaT = 100; % in msecs
totalTime = 120000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus values.
% We create here a step function of neural activity, with half-cosine ramps
%  on and off. There will be ten instances.

nInstances=10;
defaultParamsInfo.nInstances=nInstances;

for ii=1:nInstances
    stepDuration=12000; % msecs
    stepOnset=(ii-1)*stepDuration+100; % msecs
    rampDuration=3000; % msecs
    
% the square wave step
stimulusStruct.values(ii,:)=zeros(1,nTimeSamples);
stimulusStruct.values(ii,round(stepOnset/deltaT): ...
                      round(stepOnset/deltaT)+round(stepDuration/deltaT)-1)=1;
% half cosine ramp on
stimulusStruct.values(ii,round(stepOnset/deltaT): ...
                      round(stepOnset/deltaT)+round(rampDuration/deltaT)-1)= ...
                      fliplr(cos(linspace(0,pi,round(rampDuration/deltaT))/2));
% half cosine ramp off
stimulusStruct.values(ii,round(stepOnset/deltaT)+round(stepDuration/deltaT)-round(rampDuration/deltaT): ...
                      round(stepOnset/deltaT)+round(stepDuration/deltaT)-1)= ...
                      cos(linspace(0,pi,round(rampDuration/deltaT))/2);
                  
end % loop over instances building the stimulus values


%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

kernelStruct.timebase=stimulusStruct.timebase;

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct=prepareHRFKernel(kernelStruct);

%% Get the default forward model parameters
params0 = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);
params0.noiseSd = 0.05;

%% Create the packetCellArray

% We will now create a set of 10 packets. Each packet will have two
% stimulus types, each with a different expected amplitude and tau.

stimLabels={'stimA','stimB'};
stimTypes=[1 2 1 2 1 2 1 2 1 2]';
stimMeansAmplitude=[1,0.5];
nPackets=10;

for pp=1:nPackets
    
    % ensure that the packet from the pior loop is empty
    thePacket=[];
    
    % start the packet assembly
    thePacket.stimulus = stimulusStruct;
    thePacket.kernel = [];%kernelStruct;
    thePacket.metaData = [];
    
    % Randomize the order of the stimTypes
    ix = randperm(10);
    thePacket.stimulus.metaData.stimTypes=stimTypes(ix);
    thePacket.stimulus.metaData.stimLabels=stimLabels;
    
    % Create some params to define the simulated data for this packet
    paramsLocal=params0;
    for ii=1:nInstances
        paramsLocal.paramMainMatrix(ii,1)=stimMeansAmplitude(thePacket.stimulus.metaData.stimTypes(ii));
    end
    
    % Generate the simulated response
    simulatedResponseStruct = tfeHandle.computeResponse(paramsLocal,thePacket.stimulus,[],'AddNoise',true);
    
    % Add the simulated response to this packet
    thePacket.response=simulatedResponseStruct;
    
    % Add this packet to the growing cell array
    packetCellArray{pp}=thePacket;
end % loop over number of packets to be created

% Conduct the cross validation
[ trainParamsFit, trainfVals, testfVals ] = crossValidateFits( packetCellArray, tfeHandle );

