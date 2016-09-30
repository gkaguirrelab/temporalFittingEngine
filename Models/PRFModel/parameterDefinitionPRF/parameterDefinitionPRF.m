function paramStruct = parameterDefinitionPRF(nInstances,paramStructIn)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the PRF fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.

%% Unpack the params
% Here is where we would define default param values and their upper and
% lower bounds for the pRF model.
%
% MB to do some work here to allow passing of what the default params and
% their bounds should be, depending upon the stimulus properties (which
% will be specified in stimulus.metaData). Do we want this in pixels or
% degrees? Note as well that the bounds of the pRF position might need to
% go off the edge of the screen to account for voxels with receptive fields
% centered off the edge.

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