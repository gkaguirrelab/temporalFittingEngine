%% mriTemporalFittingToolboxRegistryConfig
%
% Declare the toolboxes we need for the mriTemporalFitting project and
% write them into a JSON file.  This will let us use the ToolboxToolbox to
% deliver unto us the perfect runtime environment for this project.
%
% 2016 benjamin.heasly@gmail.com

% Clear
clear;

%% Declare some toolboxes we want.
config = [ ...
    tbToolboxRecord( ...
    'type', 'git', ...
    'name', 'BrainardLabToolbox', ...
    'url', 'https://github.com/DavidBrainard/BrainardLabToolbox.git'), ...
    tbToolboxRecord( ...
    'type', 'git', ...
    'name', 'Psychtoolbox-3', ...
    'url', 'https://github.com/Psychtoolbox-3/Psychtoolbox-3.git', ...
    'subfolder','Psychtoolbox') ...   
    ];

%% Write the config to a JSON file.
configPath = 'mriTemporalFitting.json';
tbWriteConfig(config, 'configPath', configPath);

