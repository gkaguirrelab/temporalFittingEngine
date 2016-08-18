function [paramsFit,fVal,allFVals,predictedResponse] = fitResponse(obj,stimulus,responseToFit,varargin)
% [paramsFit,fVal,allFVals,predictedResponse] = fitResponse(obj,timebase,stimulus,responseToFit,varargin)
% 
% Fit method for the tmri class.  This is meant to be model independent, so
% that we only have to write it once.
%
% Inputs:
%   stimulus - cell array of stimulus descriptions
%   responseToFit - cell array of responses to fit
%
% Each entry of the cell arrays is for one run, and the fit is to find the
% parameters that minimize the average fit error, taken over the runs.
% Doing it this way allows us to more easily use this routine for cross
% validation.
% 
% Optional key/value pairs
%  'DefaultParamsInfo' - a struct passed to the defaultParams method.
%    Empty matrix is default.
%  'HRF' - a cell of structures describing the HRF to be used to go from neural to BOLD response.
%    Empty matrix is default, in which case no convolution is done
%  'paramLockMatrix' - Do parameter locking according to passed matrix.
%    This matrix has the same number of columns as the parameter vector,
%    and each row contains a 1 and a -1, which locks the two corresponding
%    parameters to each other.
%
% Outputs: 
%   paramsFit: fit parameters
%   fVal: mean value of fit error, mean taken over runs.
%   allFVals: vector of fit errors for the individual runs.
%   predictedResponse: cell array of the predicted response, one for each run

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('stimulus',@iscell);
p.addRequired('responseToFit',@iscell);
p.addParameter('DefaultParamsInfo',[],@isstruct);
p.addParameter('HRF',[],@iscell);
p.addParameter('paramLockMatrix',[],@isnumeric);
p.parse(stimulus,responseToFit,varargin{:});

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
%fFunction = @obj.fitError;
USEGLOBAL = false;
if (~USEGLOBAL)
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
    paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec, ...
       p.Results.stimulus,p.Results.responseToFit,p.Results.HRF),paramsFitVec0,[],[],p.Results.paramLockMatrix,zeros([size(p.Results.paramLockMatrix,1) 1]),vlbVec,vubVec,[],options);
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective', ...
        @(modelParamsVec)obj.fitError(modelParamsVec,p.Results.stimulus,p.Results.responseToFit,p.Results.HRF),...
        'x0',paramsFitVec0,'lb',vlbVec,'ub',vubVec,'Aeq',p.Results.paramLockMatrix,'beq',zeros([size(p.Results.paramLockMatrix,1) 1]),'options',opts);
    gs = GlobalSearch;
    paramsFitVec = run(gs,problem);
end

% Get error and predicted response for final parameters
[fVal,allFVals,predictedResponse] = obj.fitError(paramsFitVec,stimulus,responseToFit,p.Results.HRF);

% Convert fit parameters for return
paramsFit = obj.vecToParams(paramsFitVec);

end



