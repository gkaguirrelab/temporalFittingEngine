function [yPupil,gammaStimulus,ySustained,yPersistent] = createPupilTemporalModel(t,stimulus,paramsIn)

%% creatPupilTemporalModel
%
% Models the pupil temporal response as a two, temporally
% overlapping components, each controlled with an amplitude and one or
% two time-constant parameters: 
%
%  Gamma -- The passed stimulus vector is convolved witha gamma function,
%    with one parameter that defines the time constant (shape)
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
%   t - a vector of time points, in seconds
%   stimulus - a vector the describes the stimulus profile. Usually this is
%      set to have an initial value of zero and maximum value of unity.
%   paramIn - a structure of parameters that define the model. These are
%             described below.
%
% Output properties:
%
%   yPupil - the pupil response time vector, in the same time domain as t
%   gammaStimulus, ySustain, yPersistent - the decomposed components of the
%     model
%
% 07-01-2016 - gka wrote it (nc contributed the supra saturating fxn)
% 07-02-2016 - refinements to the model
%
% TO BE IMPROVED:
%
%  - better argument handling, including passing a display flag
%

%% in a more complete version, this will be set as an input flag
displayFitPlot = false;

%% vec to params
% parameters of the overall model
param.startTime = paramsIn(1); % left or right tie shift for the entire model

% parameters of the initial gamma convolution
param.gammaTau = paramsIn(2);  % time constant of the transient gamma function

% parameters of the sustained response
param.sustainedAmp = paramsIn(3); % amplitude scaling of the sustained response
param.sustainedTau = paramsIn(4); % time constant of the low-pass (exponential decay) component.

% parameters of the persistent response
param.persistentAmp = paramsIn(5); % Amplitude of the persistent filter
param.persistentT50 = paramsIn(6); % time to half-peak of the super-saturating function
param.persistentAlpha = paramsIn(7);  % time constant of the decay of the super-saturating function.

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

%% combine the elements
yPupil=sum([ySustained,yPersistent],2);

%% Plot the model
if displayFitPlot
    figure;
    subplot(1,3,1);
    hold on;
    r1 = plot(t,gammaIRF);
    r2 = plot(t,sustainedMultiplier);
    r3 = plot(t,persistentIRF);
    title('Convolution and scaling functions used in the model');
    legend([r1 r2 r3], 'gamma IRF', 'sustained Multiplier', 'persistentIRF');
    legend boxoff
    hold off;
    
    subplot(1,3,2);
    hold on;
    r1 = plot(t,gammaStimulus*(-1));
    r2 = plot(t,ySustained*(-1));
    r3 = plot(t,yPersistent*(-1));
    legend([r1 r2 r3], 'gammaStimulus', 'ySustained', 'yPersistent');
    title('Responses at sequential filter stages');
    legend boxoff
    hold off;
    
    subplot(1,3,3);
    hold on;
    r1 = plot(t,yPupil*(-1));
    legend([r1], 'yPupil');
    title('Pupil response');
    legend boxoff
    hold off;
end

gribble=1;

