function [h, meanParamValues,stdErrorParamValues] = plotParams(obj,params,stimValues,bPlot,varargin)
% plotParams(obj,params,stimulus,varargin)
%
% Generates plots of parameters for TF Block Design Color Model
%
% Inputs:
%   params   : parameter structure
%   stimulus : stimulus structure
%   bPlot    : plot or not
%
% Optional key/value pairs
%   'bFitZero' - boolean specifying whether to count '0' stimuli as stimuli
%                (default true)
% Outputs:
%   h: the figure handle
%   meanParamValues: m x n matrix of mean parameters, where m is the number
%                    of parameter types and n is the number of stimuli.
%                    Essentially stores the data to plot.

p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('stimValues',@isnumeric);
p.addRequired('bPlot',@isnumeric);
p.addParameter('bFitZero',false,@islogical);
p.parse(params,stimValues,bPlot,varargin{:});

% Optionally remove stimulus values of 0 (some analyses might want to leave them in)
if ~p.Results.bFitZero
    stimValues = stimValues(stimValues~=0);
end

% number of types of parameters, e.g. amp, tau2
numParamTypes = length(params.paramNameCell);

% get unique stimulus values
[uniqueStimValues,~] = unique(stimValues);

% for each parameter type and stimulus value
for i = 1:numParamTypes
    for j = 1:length(uniqueStimValues)
       % go into the parameter matrix, find the appropriate column, and
       % pull out all positions with a given unique stim value
       meanParamValues(i,j) = mean(params.paramMainMatrix(uniqueStimValues(j)==stimValues,i));
       stdErrorParamValues(i,j) = std(params.paramMainMatrix(uniqueStimValues(j)==stimValues,i)) ...
                                      ./length(params.paramMainMatrix(uniqueStimValues(j)==stimValues,i));
    end    
end

if bPlot == 1
    h = figure;
    set(h,'Position',[398 380 1099 664]);
    for i = 1:numParamTypes
        subplot(1,numParamTypes,i)
        plot(uniqueStimValues,meanParamValues(i,:),'-k^','MarkerSize',10,'LineWidth',1.5); axis square;
        set(gca,'Xscale','log'); set(gca,'Xtick',uniqueStimValues);
        ylabel(params.paramNameCell(i)); xlabel('Temporal frequency');
        set(gca,'FontSize',15);
    end
else
    h = [];
end

end