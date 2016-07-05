function [paramStruct,fval]= fitNeuralParams(stimMatrix,t,stimRef,data,paramStruct)
      
% set options      
options = optimoptions('fmincon','Diagnostics','on','Display','iter','Algorithm','active-set','MaxIter',100);

% grab neural parameters
prmVec0 = paramStruct.neuralParams;

f = @(prmVec)forwardModelObjectiveFunction(stimMatrix,t,stimRef,data,prmVec,paramStruct);

% set bounds
vlb = repmat(-10,[length(prmVec0) 1]);
vub = repmat(10,[length(prmVec0) 1]);

% fit parameters will also come out in a vector
[paramVec, fval] = fmincon(f,prmVec0,[],[],[],[],vlb,vub,[],options);

% we are fitting amplitudes for now
paramStruct.neuralParams = paramVec;
      
gribble = 1;