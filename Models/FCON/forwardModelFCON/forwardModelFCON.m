function [modelResponseStruct] = forwardModelFCON(params,stimulusStruct)
%% forwardModelFCON
%
% The passed parameter is a single value, termed effectiveContrast.
% The stimulusStuct contains the field "fcon" which in turn has the fields:
%
%    fcon.contrastbase - the n contrast levels at which an effective
%       contrast value is available
%    fcon.paramLookUpMatrix - a m x n matrix where m is the number of
%       parameters in the expanded description of the data, and n is the
%       number of effective contrast levels.
%    fcon.modelObjHandle - an object handle that identifies the model
%       subclass to be used to calculate the forward model using the
%       expanded parameter set.
%    fcon.defaultParams - the default params for the model subclass
%
% The contrastbase and paramLookUpMatrix is a look-up table
%  for relating an effective contrast value to an expanded set of
%  parameters.
%

%% Obtain the params
effectiveContrastVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'effectiveContrast'));

%% Define basic model features

% derive some basic properties of the stimulus values
nInstances=size(stimulusStruct.values,1);

% Set up the matrix to hold the expanded parameter set
nExpandedParameters=size(stimulusStruct.fcon.paramLookUpMatrix,1);
expandedParamMatrix=zeros(nInstances,nExpandedParameters);

% We loop through each column of the stimulus matrix and find the closest
% effective contrast in the lookup table, building a matrix of expanded
% parameters.
for ii=1:nInstances
    
    % find the closest value in fcon.contrast base to the effectiveContrast
    % for this instance
    [smallestDifference, closestIndex] = min(abs(stimulusStruct.fcon.contrastbase - effectiveContrastVec(ii)));
    
    expandedParamMatrix(ii,:)=stimulusStruct.fcon.paramLookUpMatrix(:,closestIndex);
end

% place the parameter values corresponding to the closest effective
% contrast into the default parameter structure
subclassParams=stimulusStruct.fcon.defaultParams;
subclassParams.paramMainMatrix = expandedParamMatrix;

% Obtain the model response struct for the parameter matrix
modelResponseStruct = stimulusStruct.fcon.modelObjHandle.computeResponse(subclassParams,stimulusStruct,[],'AddNoise',false);

end % function