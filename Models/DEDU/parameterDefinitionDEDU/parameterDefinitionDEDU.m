function paramStruct = parameterDefinitionDEDU(nInstances, varargin)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the DEDU fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addRequired('nInstances',@isnumeric);
p.parse(nInstances,varargin{:});

%% Unpack the params
%    These parameters are active for modeling of fMRI time series data:
%      amplitude - amplitude.
%      delay - delay in the onset of the step function in seconds
%      duration - duration of the step function in seconds
%

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    'duration',...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,2) = 1.0.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(-10000,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(0.01,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(10000,[nInstances 1]);
paramStruct.vub(:,2) = repmat(8,[nInstances 1]);

end