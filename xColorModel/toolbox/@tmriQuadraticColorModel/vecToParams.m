function params = vecToParams(obj,x)

params.Qvec(1:5) = x(1:5)';
params.crfAmp = x(6);
params.crfExponent = x(7);
params.crfSemi = x(8);
params.expFalloff = x(9);

%obj.simulateParams = params;

end