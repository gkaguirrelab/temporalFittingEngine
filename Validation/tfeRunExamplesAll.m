function status = tfeRunExamplesAll
%% Run all the examples in the tfe tree
%
% Syntax:
%     tfeRunExamplesAll
%
% Description:
%     Run all the examples in the tfe tree,
%     excepthose that contain a line of the form
%     "% ETTBSkip"
%
% Inputs:
%    None.
%
% Outputs:
%    status    - 1 if all examples run OK, 0 otherwise.
%
% Optional key/value pairs:
%    None.
%
% See also:
%   ieValidateFullAll, ieRunTutorialsAll

% History:
%   01/17/18  dhb  Wrote it.

[~, functionStatus] = ExecuteExamplesInDirectory(tbLocateToolbox('temporalFittingEngine'),'verbose',false);
status = all(functionStatus ~= -1);