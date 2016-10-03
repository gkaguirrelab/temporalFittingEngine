function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Fit method for the tfe class.  This is meant to be model independent, so
% that we only have to write it once.
%
% Inputs:
%   thePacket: a valid packet
%
% Optional key/value pairs
%  'defaultParamsInfo' - struct (default empty).  This is passed to the defaultParams method.
%  'paramLockMatrix' - matrix (default empty). If not empty, do parameter locking according
%    to passed matrix. This matrix has the same number of columns as the
%    parameter vector, and each row contains a 1 and a -1, which locks the
%    two corresponding parameters to each other.
%  'searchMethod - string (default 'fmincon').  Specify search method
%    'fmincon' - Use fmincon
%    'global' - Use global search
%    'linearRegression' - rapid estimation of simplified models with only
%                         an amplitude parameter
%
% Outputs:
%   paramsFit: fit parameters
%   fVal: mean value of fit error, mean taken over runs.
%   predictedResponse: big vector containing the fit response

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('thePacket',@isstruct);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('paramLockMatrix',[],@isnumeric);
p.addParameter('searchMethod','fmincon',@ischar);
p.addParameter('errorType','rmse',@ischar);
p.parse(thePacket,varargin{:});

% Check packet validity
if (~obj.isPacket(thePacket))
    error('The passed packet is not valid for this model');
else
    switch (obj.verbosity)
        case 'high'
            fprintf('valid\n');
    end
end

%% Set initial values and reasonable bounds on parameters
[paramsFit0,vlb,vub] = obj.defaultParams('defaultParamsInfo',p.Results.defaultParamsInfo);
paramsFitVec0 = obj.paramsToVec(paramsFit0);
vlbVec = obj.paramsToVec(vlb);
vubVec = obj.paramsToVec(vub);

%% David sez: "Fit that sucker"
switch (obj.verbosity)
    case 'high'
        fprintf('Fitting.');
end

switch (p.Results.searchMethod)
    case 'fmincon'
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
        paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec, ...
            thePacket),paramsFitVec0,[],[],p.Results.paramLockMatrix,zeros([size(p.Results.paramLockMatrix,1) 1]),vlbVec,vubVec,[],options);
    case 'global'
        % Slow but might work better
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective', ...
            @(modelParamsVec)obj.fitError(modelParamsVec,thePacket.stimulus,thePacket.response,thePacket.kernel),...
            'x0',paramsFitVec0,'lb',vlbVec,'ub',vubVec,'Aeq',p.Results.paramLockMatrix,'beq',zeros([size(p.Results.paramLockMatrix,1) 1]),'options',opts);
        gs = GlobalSearch;
        paramsFitVec = run(gs,problem);
    case 'linearRegression'
        % linear regression can be used only when the paramsFit0 has only
        % a single parameter.
        if ~length(paramsFit0.paramNameCell)==1
            error('Linear regression can only be applied in the case of a single model parameter')
        end
        %  Warn if the parameter is not called "amplitude".
        if ~(min(paramsFit0.paramNameCell{1}=='amplitude')==1)
            warning('Only amplitude parameters are suitable for linear regression')
        end
        % Take the stimulus.values as the regression matrix
        regressionMatrixStruct=thePacket.stimulus;
        % Convolve the rows of stimulus values by the kernel
        regressionMatrixStruct = obj.applyKernel(regressionMatrixStruct,thePacket.kernel,varargin{:});
        % Downsample regressionMatrixStruct to the timebase of the response
        regressionMatrixStruct = obj.resampleTimebase(regressionMatrixStruct,thePacket.response.timebase,varargin{:});
        % Perform the regression
        X=regressionMatrixStruct.values';
        y=thePacket.response.values';
        paramsFitVec=X\y;
    otherwise
        error('Do not know how to fit that sucker with specified method');
end

% Get error and predicted response for final parameters
[fVal,modelResponseStruct] = obj.fitError(paramsFitVec,thePacket,'errorType',p.Results.errorType);

switch (obj.verbosity)
    case 'high'
        fprintf('\n');
        fprintf('Fit error value: %g', fVal);
        fprintf('\n');
end

% Convert fit parameters for return
paramsFit = obj.vecToParams(paramsFitVec);

end



