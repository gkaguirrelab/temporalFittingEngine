function paramStruct = parameterDefinitionTPUP(nInstances, varargin)
% paramStruct = paramCreateBDCM(nStimuli)
%
% Create a default parameters structure for the two component
% step-function pupil model.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addRequired('nInstances',@isnumeric);
p.parse(nInstances,varargin{:});

% Parameters:
% startTime - time (in seconds) that the initial transient pupil response begins
% gammaTau - time constant of the Gamma function (msecs)
% sustainedAmp - scaling of the sustained component
% sustainedTau - time constant of the low-pass (exponential decay) component (seconds) of the sustained response
% persistentAmp - amplitude scaling of the persistent response
% persistentT50 - time to half-peak of the super-saturating function (seconds)
% persistentAlpha - time constant of the decay of the super-saturating function (seconds).




% cell for labeling each parameter column
paramStruct.paramNameCell = {...
    'delay',...
    'gammaTau', ...
    'exponentialTau', ...
    'amplitudeTransiet', ...
    'amplitudeSustained', ...
    'amplitudePersistent', ...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 100.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,2) = 150.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,3) = 100.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,4) = -100.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,5) = -100.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,6) = -100.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(0,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(100,[nInstances 1]);
paramStruct.vlb(:,3) = repmat(10,[nInstances 1]);
paramStruct.vlb(:,4) = repmat(-2000,[nInstances 1]);
paramStruct.vlb(:,5) = repmat(-2000,[nInstances 1]);
paramStruct.vlb(:,6) = repmat(-2000,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(300,[nInstances 1]);
paramStruct.vub(:,2) = repmat(300,[nInstances 1]);
paramStruct.vub(:,3) = repmat(200,[nInstances 1]);
paramStruct.vub(:,4) = repmat(0,[nInstances 1]);
paramStruct.vub(:,5) = repmat(0,[nInstances 1]);
paramStruct.vub(:,6) = repmat(0,[nInstances 1]);



end