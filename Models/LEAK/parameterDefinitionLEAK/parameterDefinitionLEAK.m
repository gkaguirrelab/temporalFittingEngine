function paramStruct = parameterDefinitionLEAK(nInstances)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the LEAK fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      timeConstant - time constant (in seconds) of the decaying exponential that
%        defines the leaky integrator
%      kappa - multiplicative adjustment of the adaptive effect of the
%        leaky integrator on the incoming signal

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    'timeConstant',...
    'kappa',...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 1.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,2) = 12.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,3) = 0.5.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(-realmax,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(.5, [nInstances 1]);
paramStruct.vlb(:,3) = repmat(0,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(realmax,[nInstances 1]);
paramStruct.vub(:,2) = repmat(60,[nInstances 1]);
paramStruct.vub(:,3) = repmat(1000,[nInstances 1]);

end