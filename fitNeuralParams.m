function [paramStruct,fval]= fitNeuralParams(stimMatrix,t,data,paramStruct)
      
% set options      
options = optimoptions('fmincon','Diagnostics','on','Display','iter','Algorithm','active-set','MaxIter',100);

% grab neural parameter matrix
prmVec0matrix = paramStruct.neuralParams;

% turn it into a vector so it can be passed into the objective function
prmVec0 = prmVec0matrix(:);

f = @(prmVec)forwardModelObjectiveFunction(stimMatrix,t,data,prmVec,paramStruct);

% make sure upper and lower bounds are also vectors
vlb = repmat([0 0 0.001],[size(prmVec0matrix,1) 1]);
vlb = vlb(:);
vub = repmat([100 1 5],[size(prmVec0matrix,1) 1]);
vub = vub(:);

% fit parameters will also come out in a vector
[paramVec, fval] = fmincon(f,prmVec0,[],[],[],[],vlb,vub,[],options);

% we are fitting amplitudes for now
paramStruct.neuralParams = reshape(paramVec,size(prmVec0matrix));
      
gribble = 1;