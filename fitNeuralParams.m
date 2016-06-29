function [paramStruct,fval]= fitNeuralParams(stimMatrix,t,data,paramStruct)
      
% set options      
options = optimoptions('fmincon','Diagnostics','on','Display','iter','Algorithm','active-set','MaxIter',100);

% convert parameter struct to initial vector
prmVec0 = paramStruct.Amplitude;

f = @(prmVec)forwardModelObjectiveFunction(stimMatrix,t,data,prmVec,paramStruct);

vlb = [];
vub = [];

[paramVec, fval] = fmincon(f,prmVec0,[],[],[],[],vlb,vub,[],options);

% we are fitting amplitudes for now
paramStruct.Amplitude = paramVec;
      
gribble = 1;