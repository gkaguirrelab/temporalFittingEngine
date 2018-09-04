function paramStruct = parameterDefinitionBTRM(nInstances, varargin)
% paramStruct = parameterDefinitionBTRM(nInstances)
%
% Create a default parameters structure for the BTRM fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.
% If passed, the paramMainMtrix (which are the initial values), vlb,
%   or vub will be used instead of the default values.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('nInstances',@isnumeric);
p.addParameter('use_dCTS',false,@islogical);
p.addParameter('initialValues',[],@isnumeric);
p.addParameter('vlb',[],@isnumeric);
p.addParameter('vub',[],@isnumeric);
p.parse(nInstances,varargin{:});

% Parameters:
% amplitude - multiplicative scaling of the stimulus.
% tauGammaIRF - time constant of the neural gamma IRF in msecs. A
%    value of 50 - 100 msecs was found in early visual areas.
%    epsilonCompression - compressive non-linearity of response.
%    Reasonable bouds are [0.1:1]. Not used if dCTS model evoked. A
%    value of 0.27 was found for area V1.
% tauInhibitoryTimeConstant - time constant (seconds) of the leaky
%    integration of neural signals that produces delayed adaptation.
% kappaInhibitionAmplitude - multiplicative scaling of the inhibitory
%    effect.
% tauExpTimeConstant - time constant of the low-pass (exponential
%    decay) component (in secs). Reasonable bounds [100:100000]
% nCompression - compressive non-linearity parameter. Reasonable
%    bounds [1:3], where 1 is no compression.
% divisiveSigma - Adjustment factor to the divisive temporal
%    normalization. Found to be ~0.1 in V1. Set to 1 to remove its effect.

% cell for labeling each parameter column
if p.Results.use_dCTS
    paramStruct.paramNameCell = { ...
        'amplitude_CTS',...
        'tauGammaIRF_CTS',...
        'tauInhibitoryTimeConstant_LEAK',...
        'kappaInhibitionAmplitude_LEAK'...
        'nCompression_dCTS',...
        'tauExpTimeConstant_dCTS',...
        'divisiveSigma_dCTS',...
        };
else
    paramStruct.paramNameCell = { ...
        'amplitude_CTS',...
        'tauGammaIRF_CTS',...
        'tauInhibitoryTimeConstant_LEAK',...
        'kappaInhibitionAmplitude_LEAK'...
        'epsilonCompression_CTS',...
        };
end

% initial values
if isempty(p.Results.initialValues)
    paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);  % amplitude_CTS
    paramStruct.paramMainMatrix(:,2) = 50.*ones([nInstances 1]);    % tauGammaIRF_CTS
    paramStruct.paramMainMatrix(:,3) = 8.*ones([nInstances 1]); % tauInhibitoryTimeConstant_LEAK
    paramStruct.paramMainMatrix(:,4) = .05.*ones([nInstances 1]); % kappaInhibitionAmplitude_LEAK
    if p.Results.use_dCTS
        paramStruct.paramMainMatrix(:,5) = 1.8.*ones([nInstances 1]);    % nCompression_dCTS
        paramStruct.paramMainMatrix(:,6) = 0.1.*ones([nInstances 1]);    % tauExpTimeConstant_dCTS
        paramStruct.paramMainMatrix(:,7) = 0.1.*ones([nInstances 1]);    % divisiveSigma_dCTS
    else
        paramStruct.paramMainMatrix(:,5) = 0.27.*ones([nInstances 1]);    % epsilonCompression_CTS
    end    
else % use passed initial values
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.paramMainMatrix(:,ii) = p.Results.initialValues(ii).*ones([nInstances 1]);
    end
end


% set lower bounds
if isempty(p.Results.vlb)
    paramStruct.vlb(:,1) = repmat(-100,[nInstances 1]);
    paramStruct.vlb(:,2) = repmat(20,[nInstances 1]);
    paramStruct.vlb(:,3) = repmat(1,[nInstances 1]);
    paramStruct.vlb(:,4) = repmat(0,[nInstances 1]);
    if p.Results.use_dCTS
        paramStruct.vlb(:,5) = repmat(1.8,[nInstances 1]);
        paramStruct.vlb(:,6) = repmat(0.1,[nInstances 1]);
        paramStruct.vlb(:,7) = repmat(0.1,[nInstances 1]);
    else
        paramStruct.vlb(:,5) = repmat(0.1,[nInstances 1]);
    end
else % used passed lower bounds
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.vlb(:,ii) = p.Results.vlb(ii).*ones([nInstances 1]);
    end
end

% set upper bounds
if isempty(p.Results.vub)
    paramStruct.vub(:,1) = repmat(100,[nInstances 1]);
    paramStruct.vub(:,2) = repmat(1000,[nInstances 1]);
    paramStruct.vub(:,3) = repmat(60,[nInstances 1]);
    paramStruct.vub(:,4) = repmat(1,[nInstances 1]);
    if p.Results.use_dCTS
        paramStruct.vub(:,5) = repmat(1.8,[nInstances 1]);
        paramStruct.vub(:,6) = repmat(0.1,[nInstances 1]);
        paramStruct.vub(:,7) = repmat(0.1,[nInstances 1]);
    else
        paramStruct.vub(:,5) = repmat(1,[nInstances 1]);
    end
else % used passed upper bounds
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.vlb(:,ii) = p.Results.vlb(ii).*ones([nInstances 1]);
    end
end

end