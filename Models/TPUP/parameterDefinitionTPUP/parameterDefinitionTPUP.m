function paramStruct = parameterDefinitionTPUP(nInstances, varargin)
% paramStruct = paramCreateBDCM(nStimuli)
%
% Create a default parameters structure for the two component
% step-function pupil model.
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
% delay - time to shift the model to the right (msecs)
% gammaTau - time constant of the Gamma function (msecs)
% exponentialTau - time constant of the persistent component (seconds)
% amplitudeTransiet - scaling of the transient component in (%change*secs)
% amplitudeSustained - scaling of the transient component in (%change*secs)
% amplitudePersistent - scaling of the transient component in (%change*secs)




% cell for labeling each parameter column
paramStruct.paramNameCell = {...
    'delay',...
    'gammaTau', ...
    'exponentialTau', ...
    'amplitudeTransiet', ...
    'amplitudeSustained', ...
    'amplitudePersistent', ...
    };

% initial values
if isempty(p.Results.initialValues)
    paramStruct.paramMainMatrix(:,1) = 200.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,2) = 200.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,3) = 10.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,4) = -10.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,5) = -25.*ones([nInstances 1]);
    paramStruct.paramMainMatrix(:,6) = -25.*ones([nInstances 1]);
else
    for ii=1:6
        paramStruct.paramMainMatrix(:,ii) = p.Results.initialValues(ii).*ones([nInstances 1]);
    end
end

% set lower bounds
if isempty(p.Results.vlb)
    paramStruct.vlb(:,1) = repmat(0,[nInstances 1]);
    paramStruct.vlb(:,2) = repmat(100,[nInstances 1]);
    paramStruct.vlb(:,3) = repmat(1,[nInstances 1]);
    paramStruct.vlb(:,4) = repmat(-2000,[nInstances 1]);
    paramStruct.vlb(:,5) = repmat(-2000,[nInstances 1]);
    paramStruct.vlb(:,6) = repmat(-2000,[nInstances 1]);
else
    for ii=1:6
        paramStruct.vlb(:,ii) = p.Results.vlb(ii).*ones([nInstances 1]);
    end
end

% set upper bounds
if isempty(p.Results.vub)
    paramStruct.vub(:,1) = repmat(500,[nInstances 1]);
    paramStruct.vub(:,2) = repmat(350,[nInstances 1]);
    paramStruct.vub(:,3) = repmat(30,[nInstances 1]);
    paramStruct.vub(:,4) = repmat(0,[nInstances 1]);
    paramStruct.vub(:,5) = repmat(0,[nInstances 1]);
    paramStruct.vub(:,6) = repmat(0,[nInstances 1]);
else
    for ii=1:6
        paramStruct.vub(:,ii) = p.Results.vub(ii).*ones([nInstances 1]);
    end
end



end