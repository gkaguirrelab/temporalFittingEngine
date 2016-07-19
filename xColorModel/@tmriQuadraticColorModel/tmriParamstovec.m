function x = tmriParamsToVec(obj)

% Take the parameter structure into a vector
for i = 1:5
    x(i) = obj.params.Qvec(i);
end
x(6) = obj.params.crfAmp;
x(7) = obj.params.crfExponent;
x(8) = obj.params.crfSemi;
x(9) = obj.params.expFalloff;

end