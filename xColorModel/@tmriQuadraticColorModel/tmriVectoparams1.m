function obj = tmriVecToParams(obj,x)

obj.params.Qvec(1:5) = x(1:5);
obj.params.crfAmp = x(6);
obj.params.crfExponent = x(7);
obj.params.crfSemi = x(8);
obj.params.expFalloff = x(9);

end