function predTC = makePredTC(stim,x0,y0,s);

[meshX , meshY] = meshgrid(1:size(stim,2),1:size(stim,1));
f = exp (-((meshY-y0).^2 + (meshX-x0).^2) ./ (2*s.^2));
for i =1:size(stim,3)
    S = stim(:,:,i).*f;
    gStim(i) = sum(S(:));
end
HRF = doubleGammaHrf(.8,[6 16],[1 1],1/6,33);
predTC = filter(HRF,1,double(gStim));
end

