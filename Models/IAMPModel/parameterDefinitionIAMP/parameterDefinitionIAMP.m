function paramStruct = parameterDefinitionIAMP(nInstances)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the IAMP fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    };

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(realmin,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(realmax,[nInstances 1]);

end