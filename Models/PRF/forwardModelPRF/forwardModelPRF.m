function [modelResponseStruct] = forwardModelPRF(params,stimulusStruct)
%% forwardModelPRF
%
% This function creates a model of neural response given a movie of
% stimulus input, a vector of time points, and a set of parameters.
%
%
% The approach is inspired by:
%

%% Obtain the params
xPosition=params.paramMainMatrix(:,strcmp(params.paramNameCell,'xPosition'));
yPosition=params.paramMainMatrix(:,strcmp(params.paramNameCell,'yPosition'));
sigmaSize=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sigmaSize'));

%% Define basic model features

% derive some basic properties of the stimulus values
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseVector=zeros(modelLength);

% Gaussian windowed sum of the stimulus at the x and y position.
[meshX , meshY] = meshgrid(1:size(stimulusStruct.values,2),1:size(stimulusStruct.values,1));
f = exp (-((meshY-xPosition).^2 + (meshX-yPosition).^2) ./ (2*sigmaSize.^2));
for i =1:size(stimulusStruct.values,3)
    S = stimulusStruct.values(:,:,i).*f;
    responseVector(i) = sum(S(:));
end

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=responseVector;

end