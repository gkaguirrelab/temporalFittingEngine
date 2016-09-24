function paramStruct = parameterDefinitionBTRM(nStimuli)
% paramStruct = paramCreateBDCM(nStimuli)
%
% Create a default parameters structure for the BTRM fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tauExponentialDecay - time constant of the low-pass (exponential
%        decay) component. Reasonable bounds [0.0001:0.1]
%
%    These parameters operate at neural timescales, so may be held fixed in
%    the modeling of fMRI data:
%      tauNeuralIRF - time constant of the neural IRF (in seconds). A
%        typical valye might be 0.005 secs (5 msecs)
%      epsilonCompression - compressive non-linearity parameter. Reasonable
%        bounds [0.1:1]

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    'tauExponentialDecay',...
    'tauNeuralIRF',...
    'epsilonCompression'...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 0.5.*ones([nStimuli 1]);
paramStruct.paramMainMatrix(:,2) = 0.001.*ones([nStimuli 1]);
paramStruct.paramMainMatrix(:,3) = 0.005.*ones([nStimuli 1]);
paramStruct.paramMainMatrix(:,4) = 0.35.*ones([nStimuli 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(-10,[nStimuli 1]);
paramStruct.vlb(:,2) = repmat(0.0001,[nStimuli 1]);
paramStruct.vlb(:,3) = repmat(0.005,[nStimuli 1]);
paramStruct.vlb(:,4) = repmat(0.35,[nStimuli 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(10,[nStimuli 1]);
paramStruct.vub(:,2) = repmat(1,[nStimuli 1]);
paramStruct.vub(:,3) = repmat(0.005,[nStimuli 1]);
paramStruct.vub(:,4) = repmat(0.35,[nStimuli 1]);

end