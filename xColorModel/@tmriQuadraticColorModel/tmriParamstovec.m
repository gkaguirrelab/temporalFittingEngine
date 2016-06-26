function x = tmriParamstovec(obj);

% Take the parameter structure into a vector
for i = 1:6
    x(i) = obj.params.Q(i);
end
x(7) = obj.params.crfAmp;
x(8) = obj.params.crfExponent;
x(9) = obj.params.crfSemi;
x(10) = obj.params.expFalloff;

end