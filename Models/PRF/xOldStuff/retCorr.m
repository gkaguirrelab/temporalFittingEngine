function minCorr = retCorr(stim,TC,paramsVec)
% minCorr = retCorr(stim,TC,paramsVec)
%   Detailed explanation goes here
x = paramsVec(1);
y = paramsVec(2);
s = paramsVec(3);
pTC = makePredTC(stim,x,y,s);
CORR = corr(TC',pTC');
minCorr = 1-CORR;
end

