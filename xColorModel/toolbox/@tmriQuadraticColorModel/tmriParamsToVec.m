function x = tmriParamsToVec(obj,varagin)

params = obj.simulateParams;

% Take the parameter structure into a vector
for i = 1:5
    x(i) = params.Qvec(i);
end
x(6) = params.crfAmp;
x(7) = params.crfExponent;
x(8) = params.crfSemi;
x(9) = params.expFalloff;

end