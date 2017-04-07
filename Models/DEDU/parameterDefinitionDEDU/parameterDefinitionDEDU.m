function paramStruct = parameterDefinitionDEDU(nInstances, varargin)
% paramStruct = paramCreateBDCM(nInstances)
%
% Create a default parameters structure for the DEDU fMRI modeling.
% This includes default parameters plus lower and upper bounds,
% as well as a field with parameter names.
%
% If passed, the paramMainMtrix (which are the initial values), vlb,
%   or vub will be used instead of the default values.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addRequired('nInstances',@isnumeric);
p.addParameter('initialValues',[],@isnumeric);
p.addParameter('vlb',[],@isnumeric);
p.addParameter('vub',[],@isnumeric);
p.parse(nInstances,varargin{:});

% Parameters:
%  amplitude - a multiplicative scaling applied to the shape of BOLD fMRI
%    response after it is normalized to have unit height.
%  duration - duration of the neural step function in seconds
%

% cell for labeling each parameter column
paramStruct.paramNameCell = { ...
    'amplitude',...
    'duration',...
    };

% initial values
if isempty(p.Results.initialValues)
    paramStruct.paramMainMatrix(:,1) = 1.0.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,2) = 1.0.*ones([nInstances 1]);
else % use passed initial values
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.paramMainMatrix(:,ii) = p.Results.initialValues(ii).*ones([nInstances 1]);
    end
end

% set lower bounds
if isempty(p.Results.vlb)
    paramStruct.vlb(:,1) = repmat(-10000,[nInstances 1]);
    paramStruct.vlb(:,2) = repmat(0.01,[nInstances 1]);
else % used passed lower bounds
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.vlb(:,ii) = p.Results.vlb(ii).*ones([nInstances 1]);
    end
end

% set upper bounds
if isempty(p.Results.vub)
    paramStruct.vub(:,1) = repmat(10000,[nInstances 1]);
    paramStruct.vub(:,2) = repmat(8,[nInstances 1]);
else % used passed upper bounds
    for ii=1:length(paramStruct.paramNameCell)
        paramStruct.vub(:,ii) = p.Results.vub(ii).*ones([nInstances 1]);
    end
end

end % function