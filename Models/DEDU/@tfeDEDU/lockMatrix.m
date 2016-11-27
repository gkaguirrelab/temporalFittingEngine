function paramLockMatrix = lockMatrix(obj,params,stimulusStruct,varargin)
% paramLockMatrix = lockMatrix(obj,params,stimulus,varargin)
%
% Return parameter locking matrix
%
% Key/value pairs
%  'LockType' - What kind of locking to do
%    'vanilla' - Default sensible locking.

%% CURRENTLY NOT IMPLEMENTED

%% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('LockType','vanilla',@ischar);
p.parse(params,varargin{:});
params = p.Results.params;

% extract temp freq values from meta data field
for i = 1:length(stimulusStruct.metaData)
   tempFreqValues(i) = stimulusStruct.metaData(i).frequency; 
end

%% Create the lock matrix
switch (p.Results.LockType)
    case 'vanilla'
        paramLockMatrix = createParamLockMatrixVanilla(unique(tempFreqValues),tempFreqValues,params.matrixCols);
    otherwise
        error('Unknown lock type passed');
end

end