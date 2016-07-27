function [paramsFit,responseFit] = fitResponse(obj,timebase,stimulus,responseToFit,varargin)
% [fitParams,fitResponse] = fitResponse(obj,timebase,stimulus,responseToFit,varargin)
% 
% Fit method for the tmri class.  This is meant to be model independent, so
% that we only have to write it once.
%
% Inputs:
%   timebase - times on which data/model predictions exist
%   stimulus - description of the stimulus
%   responseToFit - temporal response to fit
% 
% Optional key/value pairs
%  'DefaultParamsInfo' - a struct passed to the defaultParams method.
%    Empty matrix is default.
%  'HRF' - a structure describing the HRF to be used to go from neural to BOLD response.
%    Empty matrix is default, in which case no convolution is done

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('timebase',@isnumeric);
p.addRequired('stimulus')
p.addRequired('responseToFit',@isnumeric);
p.addParameter('DefaultParamsInfo',[],@isstruct);
p.addParameter('HRF',[],@isstruct);
p.parse(timebase,stimulus,responseToFit,varargin{:});

%% Set initial values and reasonable bounds on parameters
% Have a go at reasonable initial values
[paramsFit0,vlb,vub] = obj.defaultParams('DefaultParamsInfo',p.Results.DefaultParamsInfo);
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
    paramsFitVec = fmincon(@(modelParamsVec)fitFunction(modelParamsVec,obj, ...
        p.Results.timebase,p.Results.stimulus,p.Results.responseToFit,p.Results.HRF),paramsFitVec0,[],[],[],[],vlbVec,vubVec,[],options);
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective', ...
        @(modelParamsVec)fitFunction(modelParamsVec,obj,p.Results.timebase,p.Results.stimulus,p.Results.responseToFit,p.Results.HRF),...
        'x0',paramsFitVec0,'lb',vlbVec,'ub',vubVec,'options',opts);
    gs = GlobalSearch;
    paramsFitVec = run(gs,problem);
end

% Store the fit parameters
paramsFit = obj.vecToParams(paramsFitVec);
responseFit = obj.computeResponse(paramsFit,p.Results.timebase,p.Results.stimulus,'HRF',p.Results.HRF);

end

%% Error function for the fit
function f = fitFunction(paramsVec,obj,timebase,stimulus,responseToFit,HRF)

params = obj.vecToParams(paramsVec);
responsePredicted = obj.computeResponse(params,timebase,stimulus,'HRF',HRF);
f = sqrt(mean((responsePredicted-responseToFit).^2));
end
