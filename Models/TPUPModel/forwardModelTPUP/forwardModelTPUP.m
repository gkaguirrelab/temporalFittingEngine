function [modelResponseStruct] = forwardModelTPUP(params,stimulusStruct)

%% forwardModelTPUP
%
% Models the pupil temporal response as two, temporally
% overlapping components, each controlled with an amplitude and one or
% two time-constant parameters.
%
%  The passed stimulus vector is first convolved with a gamma function,
%    with one parameter that defines the time constant (shape). This
%    provides some of the "peakiness" of the initial transient seen in the
%    response data. There are then two components:
%
%  Sustain -- The convolved stimulus vector is subjected to
%     multiplicative scaling from an exponential decay function
%  Persistent -- The convolved stimulus vector is subjected to
%     convolution with a super-saturating function
%
% An additional parameter allows the entire model to shift forward or back
% in time relative to the data.
%



%% Unpack the params
%   Overall model timing
%     startTime -  left or right tie shift for the entire model
%   Parameters of the initial gamma convolution
%     gammaTau - time constant of the transient gamma function
%   Parameters of the sustained response
%     sustainedAmp - amplitude scaling of the sustained response
%     sustainedTau - time constant of the low-pass (exponential decay) component.
%   Parameters of the persistent response
%     persistentAmp - Amplitude of the persistent filter
%     persistentT50 - time to half-peak of the super-saturating function
%     persistentAlpha - time constant of the decay of the super-saturating function.

startTimeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'startTime'));
gammaTauVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'gammaTau'));
sustainedAmpVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedAmp'));
sustainedTauVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedTau'));
persistentAmpVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAmp'));
persistentT50Vec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentT50'));
persistentAlphaVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAlpha'));

% the model has parameters that are tuned for units of seconds, so
% we convert our timebase
timebaseMsecs=stimulusStruct.timebase;
timebaseSecs=timebaseMsecs/1000;

% derive some basic properties of the stimuls values
numInstances=size(stimulusStruct.values,1);
modelLength = length(timebaseSecs);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

%% We loop through each column of the stimulus matrix
for i=1:numInstances
    
    % grab the current stimulus
    stimulus=stimulusStruct.values(i,:)';
        
    %% Convolve the stimulus vector with a gamma function
    gammaIRF = timebaseSecs .* exp(-timebaseSecs/gammaTauVec(i));
    
    % scale to preserve total area after convolution
    gammaIRF=gammaIRF/sum(gammaIRF);
    
    % perform the convolution
    gammaStimulus = conv(stimulus,gammaIRF);
    gammaStimulus = gammaStimulus(1:length(timebaseSecs));
    
    %% Create the sustained component
    % Create the exponential low-pass function that defines the time-domain
    % properties of the sustain
    sustainedMultiplier=(exp(-1*sustainedTauVec(i)*timebaseSecs));
    
    % scale to preserve the max after multiplication
    sustainedMultiplier=sustainedMultiplier/max(sustainedMultiplier);
    
    % perform the multiplicative scaling
    ySustained = gammaStimulus.*sustainedMultiplier;
    
    % scale to make sure this component has unit amplitude prior to application
    % of the Amplitude parameter
    ySustained = (ySustained/max(ySustained))*sustainedAmpVec(i);    
    
    %% Create the persistent component
    % Create the super-saturating function that defines the persistent phase
    persistentIRF = createSuperSaturatingFunction(timebaseSecs,[persistentT50Vec(i),persistentAlphaVec(i)]);
    
    % scale to preserve total area after convolution
    persistentIRF=persistentIRF/sum(persistentIRF);
    
    % perform the convolution
    yPersistent = conv(gammaStimulus,persistentIRF);
    yPersistent = yPersistent(1:length(timebaseSecs));
    
    % scale the persistent component and apply the Amplitude parameter
    yPersistent = (yPersistent/max(yPersistent))*persistentAmpVec(i);
    
    %% Implement the temporal shift
    shiftAmount=find(timebaseSecs>=startTimeVec(i));
    shiftAmount=shiftAmount(1);
    gammaStimulus = circshift(gammaStimulus,[shiftAmount,0]);
    gammaStimulus(1:shiftAmount)=0;
    ySustained = circshift(ySustained,[shiftAmount,0]);
    ySustained(1:shiftAmount)=0;
    yPersistent = circshift(yPersistent,[shiftAmount,0]);
    yPersistent(1:shiftAmount)=0;
    
    %% express the model as constriction
    ySustained=ySustained*(-1);
    yPersistent=yPersistent*(-1);
    
    %% combine the elements and store
    yPupil=sum([ySustained;yPersistent],1);
    responseMatrix(i,:)=yPupil;
    
end % loop over stimulus instances

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);