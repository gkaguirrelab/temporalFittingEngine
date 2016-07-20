function fitToResponse(obj,responseToFit,varargin)
% fitToResponse(obj,responseToFit,varargin)
% 
% Compute method for the quadratic model. 

%% Set initial values and reasonable bounds on parameters
% Have a go at reasonable initial values
[modelParamsVec0,vlbVec,vubVec] = obj.defaultParams;

%% Fit that sucker
%
% I coded up the global search method, but it is very slow compared with
% fmincon alone, and fmincon seems to be fine.
USEGLOBAL = false;
if (~USEGLOBAL)
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
    modelParamsVecFit = fmincon(@(modelParamsVec)FitEllipseFunction(modelParamsVec,obj,responseToFit),modelParamsVec0,[],[],[],[],vlbVec,vubVec,[],options);
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective', ...
        @(modelParamsVec)FitEllipseFunction(modelParamsVec,obj,responseToFit),...
        'x0',modelParamsVec0,'lb',vlbVec,'ub',vubVec,'options',opts);
    gs = GlobalSearch;
    modelParamsVecFit = run(gs,problem);
end

% Store the fit parameters
obj.vecToParams(modelParamsVecFit);
obj.computeNeural;

end

%% Error function for the fit
function f = FitEllipseFunction(modelParamsVec,obj,responseToFit)

obj.vecToParams(modelParamsVec);
obj.computeNeural;
f = sqrt(mean((obj.neuralResponse-responseToFit).^2));
end
