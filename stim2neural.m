function neuralVec = stim2neural(stimMatrix,t,paramVecFit,paramStructFixed)

% function neuralVec = stim2neural(stimMatrix,paramVecFit,paramStructFixed)
%
% converts stimulus to neural model

paramIn = struct;
paramIn.afterResponseTiming = paramStructFixed.afterResponseTiming;
paramIn.tau1 = paramStructFixed.tau1;
paramIn.epsilon = paramStructFixed.epsilon;
paramIn.rectify = paramStructFixed.rectify; 

for i = 1:size(stimMatrix,1)
    paramIn.MRamplitude = paramVecFit(i,1);
    paramIn.ARampRelative = paramVecFit(i,2);
    paramIn.tau2 = paramVecFit(i,3);
    [oneStimNeural,~,~] = createNeuralTemporalModel(t,stimMatrix(i,:),0,paramIn);
    neuralMatrix(i,:) = oneStimNeural;
end

neuralVec = sum(neuralMatrix);

gribble = 1;