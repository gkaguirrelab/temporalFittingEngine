function h = plot(obj,inputStruct,varargin)
% plot(obj,timebase,response,varargin)
% 
% Plot the time-varying response 
%
% Key/value pairs
%   'NewWindow' - true/false (default true).  Create new window?
%   'Color' - vector (default [1 0 0]).  RGB color for plot.
%   'Marker' - character (default 'none'). Options include o.+x

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('responseStruct',@isstruct);
p.addParameter('NewWindow',true,@islogical);
p.addParameter('Color',[1 0 0],@isnumeric);
p.addParameter('Marker','none',@ischar);
p.addParameter('DisplayName','responseStruct.values',@ischar);
p.parse(inputStruct,varargin{:});

%% Make the figure
if (p.Results.NewWindow)
    h = figure; clf; hold on;
else
    h = gcf;
end

%% Plot
plot(inputStruct.timebase/1000,inputStruct.values, ...
    'Color',p.Results.Color,...
    'LineWidth',2,...
    'DisplayName',p.Results.DisplayName,...
    'Marker',p.Results.Marker);
xlabel('Time (seconds)')
ylabel('Response');