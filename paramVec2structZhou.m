function param = paramVec2structZhou(paramVec)

% function param = paramVec2structZhou(paramVec)
%
% converts parameter vector to struct

param.afterResponseTiming = paramVec(1);
param.MRamplitude = paramVec(2);
param.ARampRelative = paramVec(3);
param.tau1 = paramVec(4);
param.epsilon = paramVec(5);
param.tau2 = paramVec(6);

gribble = 1;