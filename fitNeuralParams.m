function [paramStruct,fval]= fitNeuralParams(stimMatrix,t,paramLockMatrix,data,paramStruct)
      
% set options      
options = optimoptions('fmincon','Diagnostics','on','Display','iter','Algorithm','active-set');

% grab neural parameters
prmVec0(:,1) = paramStruct.Amplitude;
prmVec0(:,2) = paramStruct.tau2;

f = @(prmVec)forwardModelObjectiveFunction(stimMatrix,t,data,prmVec,paramStruct);

% set bounds
vlb(:,1) = repmat(-3,[size(prmVec0,1) 1]);
vlb(:,2) = repmat(0.0001,[size(prmVec0,1) 1]);

vub(:,1) = repmat(3,[size(prmVec0,1) 1]);
vub(:,2) = repmat(0.01,[size(prmVec0,1) 1]);

% make parameter matrix into vector so can set equality constraints
prmVec0 = prmVec0(:);
vlb = vlb(:);
vub = vub(:);

Aeq = paramLockMatrix;
beq = zeros([size(paramLockMatrix,1) 1]);

% fit parameters will also come out in a vector
[paramVec, fval] = fmincon(f,prmVec0,[],[],Aeq,beq,vlb,vub,[],options);

% we are fitting amplitudes for now
paramVec = reshape(paramVec,[length(paramStruct.Amplitude) length(paramVec)./length(paramStruct.Amplitude)]);
paramStruct.Amplitude = paramVec(:,1);
paramStruct.tau2 = paramVec(:,2);
      
gribble = 1;