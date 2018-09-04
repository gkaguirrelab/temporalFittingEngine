function paramPrint(obj,params,varargin)
% print(obj,params,varargin)
%
% Print out the parameters to the command window
%
% Key/value pairs
%   'PrintType'
%     'parameters' (default)
 
%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('params',@isstruct);
p.addParameter('PrintType','parameters',@ischar);
p.parse(params,varargin{:});
params = p.Results.params;

% print the parameters
fprintf('Number of stimuli specified: %d\n',params.matrixRows);
fprintf('Number of parameter cols: %d\n',params.matrixCols);
for ii = 1:params.matrixCols
    fprintf('\tParameter %d: %s = %g\n',ii,params.paramNameCell{ii},params.paramMainMatrix(1,ii));
end
fprintf('Noise sd: %0.2f\n',params.noiseSd);


end