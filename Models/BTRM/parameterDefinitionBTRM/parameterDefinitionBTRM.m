function paramStruct = parameterDefinitionBTRM(nInstances)
% paramStruct = parameterDefinitionBTRM(nInstances)
%
% Create a default parameters structure for the BTRM fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - multiplicative scaling of the stimulus.
%      tauGammaIRF - time constant of the neural gamma IRF in msecs.
%      tauExpTimeConstant - time constant of the low-pass (exponential
%        decay) component (in mecs). Reasonable bounds [100:100000]
%      divisiveSigma - Adjustment factor to the divisive temporal
%        normalization. Found to be ~0.1 in V1. Set to unity to remove its effect. 
%      nCompression - compressive non-linearity parameter. Reasonable
%        bounds [1:3], where 1 is no compression.
%      tauInhibitoryTimeConstant - time constant of the leaky (exponential)
%        integration of neural signals that produces delayed adaptation.
%      kappaInhibitionAmplitude - multiplicative scaling of the inhibitory
%        effect.

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    'tauGammaIRF',...
    'tauExpTimeConstant',...
    'divisiveSigma',...
    'nCompression'...,...
    'tauInhibitoryTimeConstant',...
    'kappaInhibitionAmplitude',...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);  % amplitude
paramStruct.paramMainMatrix(:,2) = 1.*ones([nInstances 1]);    % tauGammaIRF
paramStruct.paramMainMatrix(:,3) = 0.001.*ones([nInstances 1]);  % tauExpTimeConstant
paramStruct.paramMainMatrix(:,4) = 0.001.*ones([nInstances 1]);    % divisiveSigma
paramStruct.paramMainMatrix(:,5) = 1.*ones([nInstances 1]);    % nCompression
paramStruct.paramMainMatrix(:,6) = 1.*ones([nInstances 1]); % tauInhibitoryTimeConstant
paramStruct.paramMainMatrix(:,7) = 0.*ones([nInstances 1]); % kappaInhibitionAmplitude

% set lower bounds
paramStruct.vlb(:,1) = repmat(-10,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(100,[nInstances 1]);
paramStruct.vlb(:,3) = repmat(0.001,[nInstances 1]);
paramStruct.vlb(:,4) = repmat(0.001,[nInstances 1]);
paramStruct.vlb(:,5) = repmat(1,[nInstances 1]);
paramStruct.vlb(:,6) = repmat(1,[nInstances 1]);
paramStruct.vlb(:,7) = repmat(0,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(10,[nInstances 1]);
paramStruct.vub(:,2) = repmat(5000,[nInstances 1]);
paramStruct.vub(:,3) = repmat(100000,[nInstances 1]);
paramStruct.vlb(:,4) = repmat(1,[nInstances 1]);
paramStruct.vub(:,5) = repmat(1,[nInstances 1]);
paramStruct.vub(:,6) = repmat(100000,[nInstances 1]);
paramStruct.vub(:,7) = repmat(1,[nInstances 1]);

end