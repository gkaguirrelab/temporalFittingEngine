function neuralVec = stim2neural(stimMatrix,t,stimRef,paramMatFit,paramStructFixed)

% function neuralVec = stim2neural(stimMatrix,paramMatFit,paramStructFixed)
%
% converts stimulus to neural model

paramIn = struct;
paramIn.afterResponseTiming = paramStructFixed.afterResponseTiming;
paramIn.tau1 = paramStructFixed.tau1;
paramIn.epsilon = paramStructFixed.epsilon;
paramIn.rectify = paramStructFixed.rectify; 
paramIn.ARampRelative = paramStructFixed.ARampRelative;
paramIn.tau2 = paramStructFixed.tau2;

for i = 1:size(stimMatrix,1)
    paramIn.MRamplitude = paramMatFit(find(paramStructFixed.stimTag==stimRef(i)),1);
%     paramIn.ARampRelative = paramMatFit(find(paramStructFixed.stimTag==stimRef(i)),2);
%     paramIn.tau2 = paramMatFit(find(paramStructFixed.stimTag==stimRef(i)),3);
    [oneStimNeural,~,~] = createNeuralTemporalModel(t,stimMatrix(i,:),0,paramIn);
    neuralMatrix(i,:) = oneStimNeural;
end

neuralVec = sum(neuralMatrix);

gribble = 1;