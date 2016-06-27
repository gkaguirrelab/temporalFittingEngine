function paramVec = paramStruct2vecZhou(param)

% function paramStruct2vecZhou(params)
%
% converts parameters for Zhou et al model to vector for fitting routine

paramVec = [param.afterResponseTiming param.MRamplitude param.ARampRelative ...
param.tau1 param.epsilon param.tau2];

gribble = 1;