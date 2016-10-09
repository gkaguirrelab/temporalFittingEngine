function paramLockMatrix = lockMatrix(obj,params,stimulus,varargin)
% paramLockMatrix = lockMatrix(obj,params,stimulus,varargin)
%
% Return parameter locking matrix.  Currently a stub.
%
% Key/value pairs
%  'LockType' - string (default 'none') What kind of locking to do?
%    'none' - Set lock matrix to empty, no lockig.

%% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('LockType','none',@ischar);
p.parse(params,varargin{:});
params = p.Results.params;

%% Create the lock matrix
switch (p.Results.LockType)
    case 'none'
        paramLockMatrix = [];
    otherwise
        error('Unknown lock type passed');
end

end