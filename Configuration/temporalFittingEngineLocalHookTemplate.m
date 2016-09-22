function temporalFittingEngineLocalHook
% temporalFittingEngine
%
% Configure things for working on the temporalFittingEngine project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "Template" from the filename,
% this will get run when you execute tbUse({'temporalFittingEngine'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% At the moment this is not used and is just a placeholder.

%% Say hello
fprintf('Running temporalFittingEngine local hook\n');

%% Put project toolbox onto path
%
% Specify project name and location
% projectName = 'mriTemporalFitting';
% projectUrl = 'https://github.com/gkaguirrelab/temporalFittingEngine.git';
% projectBaseDir = '/Users/Shared/Matlab/Analysis';
% 
% % declare the project git repo and two subfolders that we want on the path
% withTMRIToolbox = tbToolboxRecord( ...
%     'name', 'temporalFittingEngine', ...
%     'type', 'git', ...
%     'url', projectUrl, ...
%     'subfolder', 'TMRIToolbox');
% withBDCMToolbox = tbToolboxRecord( ...
%     'name', 'temporalFittingEngine', ...
%     'type', 'git', ...
%     'url', projectUrl, ...
%     'subfolder', 'BDCMToolbox');
% withQCMToolbox = tbToolboxRecord( ...
%     'name', 'temporalFittingEngine', ...
%     'type', 'git', ...
%     'url', projectUrl, ...
%     'subfolder', 'QCMToolbox');
% withPRFToolbox = tbToolboxRecord( ...
%     'name', 'temporalFittingEngine', ...
%     'type', 'git', ...
%     'url', projectUrl, ...
%     'subfolder', 'PRFToolbox');
% 
% % Obtain or update the git repo and add subfolders to the Matlab path
% config = [withTMRIToolbox withBDCMToolbox withQCMToolbox withPRFToolbox];
% tbDeployToolboxes('config', config, 'toolboxRoot', projectBaseDir, 'runLocalHooks', false);

%% Set preferences for project output
% %
% %outputBaseDir = '/Users/dhb/DropboxLab/IBIO_analysis';
% outputBaseDir = '/Volumes/Users1/DropboxLab/IBIO_analysis';
% 
% % Make base directory if it doesn't exist
% if (~exist(outputBaseDir))
%     mkdir(outputBaseDir);
% end
% 
% % This project's dir under the base dir
% theDir = fullfile(outputBaseDir,projectName);
% 
% % Set the preference
% setpref('temporalFittingEngine','outputBaseDir',theDir);
