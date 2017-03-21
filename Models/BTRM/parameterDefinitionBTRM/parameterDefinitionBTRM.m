function paramStruct = parameterDefinitionBTRM(nInstances)
% paramStruct = parameterDefinitionBTRM(nInstances)
%
% Create a default parameters structure for the BTRM fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
%      amplitude - multiplicative scaling of the stimulus.
%      tauGammaIRF - time constant of the neural gamma IRF in msecs. A
%        value of 50 - 100 msecs was found in early visual areas.
%      epsilonCompression - compressive non-linearity of response.
%        Reasonable bouds are [0.1:1]. Not used if dCTS model evoked. A
%        value of 0.27 was found for area V1.
%      tauInhibitoryTimeConstant - time constant (seconds) of the leaky
%        integration of neural signals that produces delayed adaptation.
%      kappaInhibitionAmplitude - multiplicative scaling of the inhibitory
%        effect.
%      tauExpTimeConstant - time constant of the low-pass (exponential
%        decay) component (in secs). Reasonable bounds [100:100000]
%      nCompression - compressive non-linearity parameter. Reasonable
%        bounds [1:3], where 1 is no compression.
%      divisiveSigma - Adjustment factor to the divisive temporal
%        normalization. Found to be ~0.1 in V1. Set to unity to remove its effect. 

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude_CTS',...
    'tauGammaIRF_CTS',...
    'epsilonCompression_CTS',...
    'tauInhibitoryTimeConstant_LEAK',...
    'kappaInhibitionAmplitude_LEAK'...
%     'tauExpTimeConstant_dCTS',...
%     'nCompression_dCTS',...
%     'divisiveSigma_dCTS',...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);  % amplitude_CTS
paramStruct.paramMainMatrix(:,2) = 50.*ones([nInstances 1]);    % tauGammaIRF_CTS
paramStruct.paramMainMatrix(:,3) = 0.27.*ones([nInstances 1]);    % epsilonCompression_CTS
paramStruct.paramMainMatrix(:,4) = 8.*ones([nInstances 1]); % tauInhibitoryTimeConstant_LEAK
paramStruct.paramMainMatrix(:,5) = .05.*ones([nInstances 1]); % kappaInhibitionAmplitude_LEAK
% paramStruct.paramMainMatrix(:,6) = 0.1.*ones([nInstances 1]);  % tauExpTimeConstant_dCTS
% paramStruct.paramMainMatrix(:,7) = 1.8.*ones([nInstances 1]);    % nCompression_dCTS
% paramStruct.paramMainMatrix(:,8) = 0.1.*ones([nInstances 1]);    % divisiveSigma_dCTS

% set lower bounds
paramStruct.vlb(:,1) = repmat(-100,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(50,[nInstances 1]); % 20
paramStruct.vlb(:,3) = repmat(0.27,[nInstances 1]); % 0.1
paramStruct.vlb(:,4) = repmat(1,[nInstances 1]); % 1
paramStruct.vlb(:,5) = repmat(0,[nInstances 1]); % 0
% paramStruct.vlb(:,6) = repmat(0.1,[nInstances 1]);
% paramStruct.vlb(:,7) = repmat(1,[nInstances 1]);
% paramStruct.vlb(:,8) = repmat(0.001,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(100,[nInstances 1]);
paramStruct.vub(:,2) = repmat(50,[nInstances 1]); % 1000
paramStruct.vub(:,3) = repmat(0.27,[nInstances 1]); % 1
paramStruct.vub(:,4) = repmat(50,[nInstances 1]); % 50
paramStruct.vub(:,5) = repmat(1,[nInstances 1]); % 1
% paramStruct.vub(:,6) = repmat(1,[nInstances 1]);
% paramStruct.vlb(:,7) = repmat(3,[nInstances 1]);
% paramStruct.vub(:,8) = repmat(1,[nInstances 1]);

end