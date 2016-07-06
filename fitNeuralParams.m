function [paramStruct,fval]= fitNeuralParams(stimMatrix,t,paramLockMatrix,data,paramStruct)
      
% set options      
options = optimoptions('fmincon','Diagnostics','on','Display','iter','Algorithm','active-set','MaxIter',100);

% grab neural parameters
prmVec0 = paramStruct.Amplitude;

f = @(prmVec)forwardModelObjectiveFunction(stimMatrix,t,data,prmVec,paramStruct);

% set bounds
vlb = repmat(-10,[length(prmVec0) 1]);
vub = repmat(10,[length(prmVec0) 1]);
Aeq = paramLockMatrix;
beq = zeros([size(paramLockMatrix,1) 1]);

% fit parameters will also come out in a vector
[paramVec, fval] = fmincon(f,prmVec0,[],[],Aeq,beq,vlb,vub,[],options);

% we are fitting amplitudes for now
paramStruct.Amplitude = paramVec;
      
gribble = 1;