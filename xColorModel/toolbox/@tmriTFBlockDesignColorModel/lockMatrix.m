function paramLockMatrix = lockMatrix(obj,params,stimulus,varargin)
% paramLockMatrix = lockMatrix(obj,params,stimulus,varargin)
%
% Return parameter locking matrix
%
% Key/value pairs
%  'LockType' - What kind of locking to do
%    'vanilla' - Default sensible locking.

%% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('LockType','vanilla',@ischar);
p.parse(params,varargin{:});
params = p.Results.params;

%% Create the lock matrix
switch (p.Results.LockType)
    case 'vanilla'
        paramLockMatrix = createParamLockMatrixVanilla(stimulus.uniqueTemporalFreq,stimulus.stimValues,params.matrixCols);
    otherwise
        error('Unknown lock type passed');
end

end