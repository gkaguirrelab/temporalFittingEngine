function [paramsFit,responseFit] = fitResponse(obj,timebase,stimulus,responseToFit,varargin)
% [fitParams,fitResponse] = fitResponse(obj,timebase,stimulus,responseToFit,varargin)
% 
% Compute method for the quadratic model. 

%% Set initial values and reasonable bounds on parameters
% Have a go at reasonable initial values
[paramsFit0,vlb,vub] = obj.defaultParams;
paramsFitVec0 = obj.paramsToVec(paramsFit0);
vlbVec = obj.paramsToVec(vlb);
vubVec = obj.paramsToVec(vub);

%% Fit that sucker
%
% I coded up the global search method, but it is very slow compared with
% fmincon alone, and fmincon seems to be fine.
USEGLOBAL = false;
if (~USEGLOBAL)
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
    paramsFitVec = fmincon(@(modelParamsVec)FitEllipseFunction(modelParamsVec,obj,timebase,stimulus,responseToFit),paramsFitVec0,[],[],[],[],vlbVec,vubVec,[],options);
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective', ...
        @(modelParamsVec)FitEllipseFunction(modelParamsVec,obj,timebase,stimulus,responseToFit),...
        'x0',paramsFitVec0,'lb',vlbVec,'ub',vubVec,'options',opts);
    gs = GlobalSearch;
    paramsFitVec = run(gs,problem);
end

% Store the fit parameters
paramsFit = obj.vecToParams(paramsFitVec);
responseFit = obj.computeResponse(paramsFit,timebase,stimulus);

end

%% Error function for the fit
function f = FitEllipseFunction(paramsVec,obj,timebase,stimulus,responseToFit)

params = obj.vecToParams(paramsVec);
responsePredicted = obj.computeResponse(params,timebase,stimulus);
f = sqrt(mean((responsePredicted-responseToFit).^2));
end
