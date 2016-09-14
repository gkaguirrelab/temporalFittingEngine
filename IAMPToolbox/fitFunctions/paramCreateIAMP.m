function paramStruct = paramCreateIAMP(nStimuli)
% paramStruct = paramCreateIAMP(nStimuli)
%
% Create a default parameters structure for the instance amplitude modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.


% cell for labeling each parameter column
paramStruct.paramNameCell = {'Amplitude'};

% initial values
paramStruct.paramMainMatrix = [];
paramStruct.paramMainMatrix(:,1) = ones([nStimuli 1]);

% set lower bounds
paramStruct.vlb(:,1) = repmat(-10,[nStimuli 1]);

% set upper bounds
paramStruct.vub(:,1) = repmat(10,[nStimuli 1]);

end