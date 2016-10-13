function paramStruct = parameterDefinitionFCON(nInstances)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the FCON fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'effectiveContrast',...
    };

% initial values
paramStruct.paramMainMatrix(:,1) = 0.*ones([nInstances 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(-1,[nInstances 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(1,[nInstances 1]);

end