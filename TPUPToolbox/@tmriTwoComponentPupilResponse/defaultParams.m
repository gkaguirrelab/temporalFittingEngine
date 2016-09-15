function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct form, use paramsToVec on the structs to
% get vector form.
%
% % Optional key/value pairs
%  'DefaultParamsInfo' - A struct passed to the defaultParams method.  This
%  struct should have a field called nStimuli, which is used by this
%  routine.  If this struct is not passed, nStimuli is set to 10.

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('DefaultParamsInfo',[],@isstruct);
p.parse(varargin{:});

%% Handle default number of stimuli
if (isempty(p.Results.DefaultParamsInfo))
    nInstances = 10;
else
    nInstances = p.Results.DefaultParamsInfo.nInstances;
end

% Use general routine that is also called by the non-object oriented
% version of this model.  For that reason, we'll need to do a little
% reformating.
paramStruct = parameterDefinitionTPUP(nInstances);
params.paramNameCell = paramStruct.paramNameCell;
params.paramMainMatrix = paramStruct.paramMainMatrix;
params.matrixRows = size(params.paramMainMatrix,1);
params.matrixCols = size(params.paramMainMatrix,2);

% Upper and lower bounds
paramsLb.paramNameCell = paramStruct.paramNameCell;
paramsLb.paramMainMatrix = paramStruct.vlb;
paramsLb.matrixRows = size(paramsLb.paramMainMatrix,1);
paramsLb.matrixCols = size(paramsLb.paramMainMatrix,2);
paramsUb.paramNameCell = paramStruct.paramNameCell;
paramsUb.paramMainMatrix = paramStruct.vub;
paramsUb.matrixRows = size(paramsUb.paramMainMatrix,1);
paramsUb.matrixCols = size(paramsUb.paramMainMatrix,2);

% Noise parameter for simulation
params.noiseSd = 0;
paramsLb.noiseSd = 0;
paramsUb.noiseSd = 0;


end