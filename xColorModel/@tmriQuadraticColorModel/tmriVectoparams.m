function obj = tmriVectoparams(obj,x)

obj.params.Q(1:6) = x(1:6);
obj.params.crfAmp = x(7);
obj.params.crfExponent = x(8);
obj.params.crfSemi = x(9);
obj.params.expFalloff = x(10);

end