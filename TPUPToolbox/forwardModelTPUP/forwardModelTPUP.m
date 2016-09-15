function [responseMatrix] = forwardModelTPUP(t,stimulusMatrix, startTimeVec, gammaTauVec, sustainedAmpVec, sustainedTauVec, persistentAmpVec, persistentT50Vec, persistentAlphaVec)

%% creatPupilTemporalModel
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
% Input properties:
%
%   t - a vector of time points, in milliseconds
%   stimulus - a matrix of stimulus instances [instances, time].
%      Usually, each instance is a vector with an initial value of zero,
%      and maximum value of unity, although this is not necessary.
%   A set of parameters that define the model. These are
%             described below.
%
% Output properties:
%
%   yPupil - the pupil response time vector, in the same time domain as t
%
% 07-01-2016 - gka wrote it (nc contributed the supra saturating fxn)
% 07-02-2016 - refinements to the model
% 09-13-2016 - converted to be used within the temporalFittingEngine
%

% derive some basic properties of the stimuls matrix and model
numInstances=size(stimMatrix,1);
modelLength = length(t);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

%% We loop through each column of the stimulus matrix
for i=1:numInstances
    
    % grab the current stimulus
    stimulus=stimulusMatrix(i,:);
    
    % grab the params for this stimulus instance
    
    % parameters of the overall model
    param.startTime = startTimeVec(i); % left or right tie shift for the entire model
    
    % parameters of the initial gamma convolution
    param.gammaTau = gammaTauVec(i);  % time constant of the transient gamma function
    
    % parameters of the sustained response
    param.sustainedAmp = sustainedAmpVec(i); % amplitude scaling of the sustained response
    param.sustainedTau = sustainedTauVec(i); % time constant of the low-pass (exponential decay) component.
    
    % parameters of the persistent response
    param.persistentAmp = persistentAmpVec(i); % Amplitude of the persistent filter
    param.persistentT50 = persistentT50Vec(i); % time to half-peak of the super-saturating function
    param.persistentAlpha = persistentAlphaVec(i);  % time constant of the decay of the super-saturating function.
    
    %% Convolve the stimulus vector with a gamma function
    gammaIRF = t .* exp(-t/param.gammaTau);
    
    % scale to preserve total area after convolution
    gammaIRF=gammaIRF/sum(gammaIRF);
    
    % perform the convolution
    gammaStimulus = conv(stimulus,gammaIRF);
    gammaStimulus = gammaStimulus(1:length(t));
    
    %% Create the sustained component
    % Create the exponential low-pass function that defines the time-domain
    % properties of the sustain
    sustainedMultiplier=(exp(-1*param.sustainedTau*t));
    
    % scale to preserve the max after multiplication
    sustainedMultiplier=sustainedMultiplier/max(sustainedMultiplier);
    
    % perform the multiplicative scaling
    ySustained = gammaStimulus.*sustainedMultiplier;
    
    % scale to make sure this component has unit amplitude prior to application
    % of the Amplitude parameter
    ySustained = (ySustained/max(ySustained))*param.sustainedAmp;    
    
    %% Create the persistent component
    % Create the super-saturating function that defines the persistent phase
    persistentIRF = createSuperSaturatingFunction(t,[param.persistentT50,param.persistentAlpha]);
    
    % scale to preserve total area after convolution
    persistentIRF=persistentIRF/sum(persistentIRF);
    
    % perform the convolution
    yPersistent = conv(gammaStimulus,persistentIRF);
    yPersistent = yPersistent(1:length(t));
    
    % scale to make sure this component has unit amplitude prior to application
    % of the Amplitude parameter
    yPersistent = (yPersistent/max(yPersistent))*param.persistentAmp;
    
    %% Implement the temporal shift
    shiftAmount=find(t>=param.startTime);
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
    yPupil=sum([ySustained,yPersistent],2);
    responseMatrix(i,:)=yPupil;
    
end % loop over stimulus instances
