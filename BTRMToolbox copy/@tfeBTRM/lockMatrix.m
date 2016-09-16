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
params = p.results.params;

% extract temp freq values from meta data field
for i = 1:length(stimulus.metaData)
   tempFreqValues(i) = stimulus.metaData(i).frequency; 
end

%% Create the lock matrix
switch (p.results.LockType)
    case 'vanilla'
        paramLockMatrix = createVanillaLockMatrixBTRM(unique(tempFreqValues),tempFreqValues,params.matrixCols);
    otherwise
        error('Unknown lock type passed');
end

end