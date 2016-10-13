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

%% Specify project-specific preferences
%
% This currently include UnitTestToolbox/RemoteDataToolbox setup
p = struct(...
    'projectName',           'temporalFittingEngine', ...                                                                     % The project's name (also the preferences group name)
    'validationRootDir',     tfeValidationDir, ...                                                                            % Directory location where the 'scripts' subdirectory resides.
    'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
    'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
    'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
    'remoteDataToolboxConfig', 'temporalFittingEngine', ...                                                                   % Struct, file path, or project name with Remote Data Toolbox configuration.
    'clonedWikiLocation',    '', ...                                                                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
    'clonedGhPagesLocation', '', ...                                                                                          % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
    'githubRepoURL',         '', ...                                                                                          % Github URL for the project. This is only used for publishing tutorials.
    'generateGroundTruthDataIfNotFound',false,...                                                                            % Flag indicating whether to generate ground truth if one is not found
    'listingScript',         'tfeValidateListAllValidationDirs', ...                                                          % Script that lists dirs to find validation scripts in
    'coreListingScript',     '', ...                                                                                          % Not used in this project
    'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
    );

generatePreferenceGroup(p);

%% Only do this if UnitTestToolbox is on the path.
%
% It might not be if we're working on a content analysis project that uses the
% temporalFittingEngine.
if exist('UnitTest','file')
    UnitTest.usePreferencesForProject(p.projectName);
end

end

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

function generatePreferenceGroup(p)
    % Remove any existing preferences for this project
    if ispref(p.projectName)
        rmpref(p.projectName);
    end
    
    % Renerate and save the project-specific preferences
    setpref(p.projectName, 'projectSpecificPreferences', p);
    fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end
