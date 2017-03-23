function paramStruct = parameterDefinitionPRF(nInstances,varargin)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the pRF fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addRequired('nInstances',@isnumeric);
p.parse(nInstances,varargin{:});

%% Unpack the params
% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'xPosition',...
    'yPosition',...
    'sigmaSize'...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 540.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,2) = 540.*ones([nInstances 1]);
paramStruct.paramMainMatrix(:,3) = 10.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(0,[nInstances 1]);
paramStruct.vlb(:,2) = repmat(0,[nInstances 1]);
paramStruct.vlb(:,3) = repmat(1,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(1080,[nInstances 1]);
paramStruct.vub(:,2) = repmat(1080,[nInstances 1]);
paramStruct.vub(:,3) = repmat(300,[nInstances 1]);

end