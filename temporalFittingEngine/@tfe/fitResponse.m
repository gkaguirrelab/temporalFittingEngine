function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Syntax:
%  [paramsFit,fVal,modelResponseStruct] = obj.fitResponse(thePacket)
%
% Description:
%   Fit method for the tfe class.  This is meant to be model independent,
%   so that we only have to write it once.
%
% Inputs:
%   thePacket             - Structure. A valid packet
%
% Optional key/value pairs
%  'defaultParamsInfo'    - Struct (default empty).  This is passed to the
%                           defaultParams method. This can be used to do
%                           custom things for the search routine. [We might
%                           try to make this go away.]
%  'defaultParams'        - Struct (default empty). Params values for
%                           defaultParams to return. In turn determines
%                           starting value for search, unless overridden by
%                           a non-empty initialParams key/value pair. [We
%                           might try to make this go away, and use the
%                           initialParams key/value pair below, which is
%                           simpler to read and think about.]
%  'intialParams'         - Struct (default empty). Params structure
%                           containing initial parameters for search. If empty,
%                           what's returned by the defaultParams method is
%                           used. This key/value pair overrides the
%                           defaultParams key/value pair, which in turn
%                           overrides what is returned by the
%                           defaultParams method.
%  'vlbParams'            - Struct (default empty). Params structure
%                           containing lower bounds for search. If empty,
%                           what's returned by the defaultParams method is
%                           used.
%  'vubParams'            - Struct (default empty). Params structure
%                           containing upper bounds for search. If empty,
%                           what's returned by the defaultParams method is
%                           used.
%  'A'                     - Matrix (default empty). Linear inequality
%                           constraint matrix.
%  'b'                    - Vector (default empty). Linear inequality
%                           constraint vector.
%  'Aeq'                  - Matrix (default empty). Linear equality
%                           constraint matrix.
%  'beq'                  - Vector (default empty). Linear equality
%                           constraint vector.
%  'nlcon'                - Handle to nonlinear constraint function
%                           (default empty)
%  'searchMethod          - String (default 'fmincon').  Specify search
%                           method:
%                              'fmincon' - Use fmincon
%                              'linearRegression' - Rapid estimation of
%                                   simplified models with only an
%                                   amplitude parameter. Parameter bounds
%                                   and constraints do not apply. Nor does
%                                   errorType or errorWeightVector.
%  'DiffMinChange'        - Double (default empty). If not empty, changes
%                           the default value of this in the fmincon
%                           optionset.
%  'fminconAlgorithm'     - String (default 'active-set'). If set to a string,
%                           passed on as algorithm in options to fmincon.
%                           Can be empty or any algorithm string understood
%                           by fmincon.
%                              [] - Use fmincon's current default algorithm
%                              'active-set' - Active set algorithm
%                              'interior-point' - Interior point algorithm.
%  'errorType'            - String (default 'rmse'). Determines what error
%                           is minimized, passed along as an option to the
%                           fitError method.
%  'errorWeightVector'    - Vector of weights to use on error for each
%                           response value. Only valid if error type is
%                           'rmse'. Passed along as an option to the
%                           fitError method.
%  'fitErrorScalar'       - Computed fit error is multiplied by this before
%                           return.  Sometimes getting the objective
%                           function onto the right scale makes all the
%                           difference in fitting. Passed along as an
%                           option to the fitError method.
%  'maxIter'              - If not empty, use this as the maximum number of
%                           interations in the search.
%  'maxFunEval'           - If not empty, use this as the maximum number of
%                           function evaluations in the search.
%
% Outputs:
%   paramsFit             - Structure. Fit parameters.
%   fVal                  - Scalar. Fit error
%   predictedResponse     - Structure. Response predicted from fit.

% History:
%  11/26/18  dhb       Added comments about key/value pairs that were not
%                      previously commented.
%  12/09/18  dhb       Comment improvements.
%  12/13/18  gka       Notes placed in git commit comments.
%  01/01/18  dhb       Allows explicit passing of initial parameters,
%                      bounds, and search constraints.
%            dhb       Code cleaning. Remove comment about 'global'
%                      searchMethod value, because it doesn't exist in the code.
%            dhb       Final computation of error was not using
%                      errorWeightVector.  Now passed on.
%  04/03/21  dhb       Allow passing of nlcon handle to fmincon with
%                      key/value pair.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('thePacket',@isstruct);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('initialParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('vlbParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('vubParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('A',[],@isnumeric);
p.addParameter('b',[],@isnumeric);
p.addParameter('Aeq',[],@isnumeric);
p.addParameter('beq',[],@isnumeric);
p.addParameter('nlcon',[]);
p.addParameter('maxIter',[],@isnumeric);
p.addParameter('maxFunEval',[],@isnumeric);
p.addParameter('searchMethod','fmincon',@ischar);
p.addParameter('DiffMinChange',[],@isnumeric);
p.addParameter('fminconAlgorithm','active-set',@(x) (isempty(x) | ischar(x)));
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
% Uses defaultParams method.
[initialParams,vlbParams,vubParams] = obj.defaultParams('defaultParamsInfo',p.Results.defaultParamsInfo,'defaultParams',p.Results.defaultParams,varargin{:});

%% Optional override of initial params and bounds
if (~isempty(p.Results.initialParams))
    initialParams = p.Results.initialParams;
end
if (~isempty(p.Results.vlbParams))
    vlbParams = p.Results.vlbParams;
end
if (~isempty(p.Results.vubParams))
    vubParams = p.Results.vubParams;
end

%% Convert initial params and bounds to vector form
initialParamsVec = obj.paramsToVec(initialParams);
vlbVec = obj.paramsToVec(vlbParams);
vubVec = obj.paramsToVec(vubParams);

%% David sez: "Fit that sucker"
switch (obj.verbosity)
    case 'high'
        fprintf('Fitting.');
end

switch (p.Results.searchMethod)
    case 'fmincon'
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off');
        if (~isempty(p.Results.fminconAlgorithm))
            options = optimset(options,'Algorithm',p.Results.fminconAlgorithm);
        end
        if ~isempty(p.Results.DiffMinChange)
            options = optimset(options,'DiffMinChange',p.Results.DiffMinChange);
        end
        if (~isempty(p.Results.maxIter))
            options = optimset(options,'MaxIter',p.Results.maxIter);
        end
        if (~isempty(p.Results.maxFunEval))
            options = optimset(options,'MaxFunEval',p.Results.maxFunEval);
        end

        % Uncomment to call a function so that parameter values are printed
        % out on each displayed iteration.  Only for deep debugging.
        % You'll probably want to set 'Display' to 'iter' above when using
        % this.
        %
        % options = optimset('OutputFcn',@outputFunction);

        % Do the actual fit
        paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec,thePacket,varargin{:}), ...
            initialParamsVec,p.Results.A,p.Results.b,p.Results.Aeq,p.Results.beq,vlbVec,vubVec,p.Results.nlcon,options);

    case 'linearRegression'
        % Linear regression can be used only when the paramsFit0 has only a
        % single parameter.  In addition, constraints,bounds,errorType and
        % errorWeightVector do not apply to this method.
        if length(initialParams.paramNameCell)~=1
            error('fitResponse:invalidLinearRegression','Linear regression can only be applied in the case of a single model parameter')
        end

        %  Warn if the parameter is not called "amplitude".
        if ~(min(initialParams.paramNameCell{1}=='amplitude')==1)
            warning('fitResponse:invalidLinearRegression','Only amplitude parameters are suitable for linear regression')
        end

        % Take the stimulus.values as the regression matrix
        regressionMatrixStruct=thePacket.stimulus;

        % Convolve the rows of stimulus values by the kernel
        regressionMatrixStruct = obj.applyKernel(regressionMatrixStruct,thePacket.kernel,varargin{:});

        % Downsample regressionMatrixStruct to the timebase of the response
        regressionMatrixStruct = obj.resampleTimebase(regressionMatrixStruct,thePacket.response.timebase,varargin{:});

        % Pull out the response from the packet
        y=thePacket.response.values';

        % Assign the regression matrix to X
        X=regressionMatrixStruct.values';

        % Detect if nans are present in the response. If so, remove these
        % timepoints from the response and corresponding locations in the
        % regression matrix
        if any(isnan(y))
            validIdx = ~isnan(y);
            y = y(validIdx);
            X = X(validIdx,:);
        end

        % Perform the regression
        paramsFitVec=X\y;
    case 'bads'
        % Do the actual fit
%         paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec,thePacket,varargin{:}), ...
%             initialParamsVec,p.Results.A,p.Results.b,p.Results.Aeq,p.Results.beq,vlbVec,vubVec,p.Results.nlcon,options);
    [paramsFitVec,fval] = bads(@(modelParamsVec)obj.fitError(modelParamsVec,thePacket,varargin{:}),...
                     initialParamsVec',vlbVec',vubVec',[],[])
    otherwise
        error('fitResponse:invalidSearchMethod','Do not know how to fit that sucker with specified method');
end

% Get error and predicted response for final parameters
[fVal,modelResponseStruct] = obj.fitError(paramsFitVec,thePacket,varargin{:});

switch (obj.verbosity)
    case 'high'
        fprintf('\n');
        fprintf('Fit error value: %g', fVal);
        fprintf('\n');
end

% Convert fit parameters for return
paramsFit = obj.vecToParams(paramsFitVec);

end

%% You can use this code to get the parameter values printed on each display iteration.
function stop = outputFunction(x, optimValues, state)

if (strcmp(state,'iter'))
    x'
end
stop = false;

end


