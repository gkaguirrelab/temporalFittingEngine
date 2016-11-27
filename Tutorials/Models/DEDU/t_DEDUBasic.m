function [  ] = t_DEDUBasic(varargin)
% function [  ] = t_DEDUBasic(varargin)
%
% Demonstrate the DElay and DUration model
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
tfeHandle = tfeDEDU('verbosity','none');

%% Temporal domain of the stimulus
deltaT = 10; % in msecs
totalTime = 200000; % in msecs
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);

%% Specify the stimulus values.
% We create here a step function of neural activity, with half-cosine ramps
%  on and off. There will be ten instances.
nInstances=10;
stepDuration=3000; % msecs
onsetDelay=100; % msecs
rampDuration=250; % msecs
spacing=16000;

defaultParamsInfo.nInstances=nInstances;

for ii=1:nInstances

    % calculate the time at which the stimulus begins
    stepOnset=(ii-1)*spacing+100; % msecs
    
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
params0.noiseSd = 0.005;

%% Create the packetCellArray

% We will now create a set of nPackets. Each packet will have two
% stimulus types, each with a different expected amplitude and tau.
stimLabels={'stimA','stimB'};
stimTypes=[1 2 1 2 1 2 1 2 1 2]';
stimMeansAmplitude=[1,0.9375];
stimMeansDelay=[0,.3];
stimMeansDuration=[3,3.5];

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
        paramsLocal.paramMainMatrix(ii,2)=stimMeansDelay(thePacket.stimulus.metaData.stimTypes(ii));
        paramsLocal.paramMainMatrix(ii,3)=stimMeansDuration(thePacket.stimulus.metaData.stimTypes(ii));
    end
    
    % Generate the simulated response
    simulatedResponseStruct = tfeHandle.computeResponse(paramsLocal,thePacket.stimulus,thePacket.kernel,'AddNoise',true);
    
    % Add the simulated response to this packet
    thePacket.response=simulatedResponseStruct;
    
tfeHandle.plot(simulatedResponseStruct,'DisplayName','Simulated');

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

%% Test the fitter
[paramsFit,fVal,modelResponseStruct] = ...
            tfeHandle.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix, ...
            'DiffMinChange',0.01);
        
%% Report the output
fprintf('Model parameter from fits:\n');
tfeHandle.paramPrint(paramsFit);
fprintf('\n');

tfeHandle.plot(modelResponseStruct,'Color',[0 1 0],'NewWindow',false,'DisplayName','model fit');
legend('show');legend('boxoff');
hold off;

end % function
