function [modelResponseStruct] = forwardModelDEDU(obj,params,stimulusStruct, hrfKernelStruct)
%% forwardModelDEDU
%
% This function creates a model of neural response given a vector of
% stimulus input, a vector of time points, and a set of parameters.
%
% The model of neural activity is a step function, controlled by three
% parameters: amplitude, delay, and duration (secs). A kernel is built
% based upon the parameters, and then the stimulus vector is convolved by
% this kernel.
%

amplitudeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude'));
durationVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'duration'));

%% Define basic model features

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
modelLength = length(stimulusStruct.timebase);

% pre-allocate the responseMatrix variable here for speed
responseMatrix=zeros(numInstances,modelLength);

% determine the deltaT of the stimulus and make sure it is regularly
% sampled (in msecs)
check = diff(stimulusStruct.timebase);
deltaT = check(1);
if (any(abs(check - check(1)) > 1e-6))
    error('Response structure timebase is not regularly sampled');
end

%% We loop through each column of the stimulus matrix
for i=1:numInstances
    
    % grab the current stimulus
    subStimulusStruct=stimulusStruct;
    subStimulusStruct.values=subStimulusStruct.values(i,:);

    % Build the kernel
    stepKernelStruct.timebase=0:deltaT:max([ceil(durationVec(i)*1000-deltaT) deltaT]);
    stepKernelStruct.values=stepKernelStruct.timebase*0+1;
    
    % Convolve the stimulus by the step kernel
    yNeuralStruct=obj.applyKernel(subStimulusStruct,stepKernelStruct);

    % Convolve the stimulus by the hrf kernel
    yNeuralStruct=obj.applyKernel(yNeuralStruct,hrfKernelStruct);

    % make the response have unit amplitude
    yNeuralStruct.values=yNeuralStruct.values./max(yNeuralStruct.values);

    % Scale by the amplitude parameter
    yNeuralStruct.values=yNeuralStruct.values.*amplitudeVec(i);
    
    %% Place yNeural into the growing neuralMatrix
    responseMatrix(i,:)=yNeuralStruct.values';
    
end % loop over rows of the stimulus matrix

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end