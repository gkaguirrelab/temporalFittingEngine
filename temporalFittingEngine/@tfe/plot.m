function h = plot(obj,responseStruct,varargin)
% plot(obj,timebase,response,varargin)
% 
% Plot the time-varying response 
%
% Key/value pairs
%   'NewWindow' - true/false (default true).  Create new window?
%   'Color' - vector (default [1 0 0]).  RGB color for plot.

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('responseStruct',@isstruct);
p.addParameter('NewWindow',true,@islogical);
p.addParameter('Color',[1 0 0],@isnumeric);
p.addParameter('DisplayName','responseStruct.values',@ischar);
p.parse(responseStruct,varargin{:});

%% Make the figure
if (p.Results.NewWindow)
    h = figure; clf; hold on;
else
    h = gcf;
end

%% Plot
plot(responseStruct.timebase,responseStruct.values, ...
    'Color',p.Results.Color,...
    'LineWidth',2,...
    'DisplayName',p.Results.DisplayName);
xlabel('Time (seconds)')
ylabel('Response');