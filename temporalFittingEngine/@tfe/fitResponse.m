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
%  'fminconAlgorithm'     - String (default 'active-set'). Passed on as
%                           algorithm in options to fmincon.
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
%           dhb        Final computation of error was not using
%                      errorWeightVector.  Now passed on.

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
p.addParameter('searchMethod','fmincon',@ischar);
p.addParameter('DiffMinChange',[],@isnumeric);
p.addParameter('fminconAlgorithm','active-set',@ischar);
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
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm',p.Results.fminconAlgorithm);
        if ~isempty(p.Results.DiffMinChange)
            options = optimset(options,'DiffMinChange',p.Results.DiffMinChange);
        end
        paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec,thePacket,varargin{:}), ...
            initialParamsVec,p.Results.A,p.Results.b,p.Results.Aeq,p.Results.beq,vlbVec,vubVec,[],options);
    
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
        
        % Perform the regression
        X=regressionMatrixStruct.values';
        y=thePacket.response.values';
        paramsFitVec=X\y;
        
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



