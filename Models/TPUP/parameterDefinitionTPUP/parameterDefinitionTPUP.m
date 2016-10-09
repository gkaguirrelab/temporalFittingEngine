function paramStruct = parameterDefinitionTPUP(nInstances)
% paramStruct = paramCreateBDCM(nStimuli)
%
% Create a default parameters structure for the two component
% step-function pupil model.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

% Parameters:
% startTime - time (in seconds) that the initial transient pupil response begins
% gammaTau - time constant of the Gamma function (msecs)
% sustainedAmp - scaling of the sustained component
% sustainedTau - time constant of the low-pass (exponential decay) component (seconds) of the sustained response
% persistentAmp - amplitude scaling of the persistent response
% persistentT50 - time to half-peak of the super-saturating function (seconds)
% persistentAlpha - time constant of the decay of the super-saturating function (seconds).

% cell for labeling each parameter column
paramStruct.paramNameCell = {'startTime','gammaTau', ...
                             'sustainedAmp', 'sustainedTau', ...
                             'persistentAmp','persistentT50','persistentAlpha'};

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 0.1.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,2) = 0.2.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,3) = 0.2.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,4) = 0.5.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,5) = 0.2.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,6) = 0.25.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,7) = 1.5.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(0.1,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(0.1,[nInstances 1]);
paramStruct.vlb(:,3) = repmat(eps,[nInstances 1]);
paramStruct.vlb(:,4) = repmat(eps,[nInstances 1]);
paramStruct.vlb(:,5) = repmat(eps,[nInstances 1]);
paramStruct.vlb(:,6) = repmat(eps,[nInstances 1]);
paramStruct.vlb(:,7) = repmat(0.5,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(0.4,[nInstances 1]);
paramStruct.vub(:,2) = repmat(0.3,[nInstances 1]);
paramStruct.vub(:,3) = repmat(1,[nInstances 1]);
paramStruct.vub(:,4) = repmat(2,[nInstances 1]);
paramStruct.vub(:,5) = repmat(1,[nInstances 1]);
paramStruct.vub(:,6) = repmat(2,[nInstances 1]);
paramStruct.vub(:,7) = repmat(6.0,[nInstances 1]);



end