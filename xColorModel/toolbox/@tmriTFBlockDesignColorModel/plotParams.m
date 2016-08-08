function [h, meanParamValues] = plotParams(obj,params,stimulus,varargin)
% plotParams(obj,param,stimulus,varargin)
%
% Generates plots of parameters for TF Block Design Color Model
%
% Inputs
%   params   : parameter structure
%   stimulus : stimulus structure
% Optional Inputs
%   bFitZero : boolean specifying whether to count '0' stimuli as stimuli
% Outputs
% h              : the figure
% meanParamValues: m x n matrix of mean parameters, where m is the number
%                  of parameter types and n is the number of stimuli.
%                  Essentially stores the data to plot.

p = inputParser;
p.addRequired('param',@isstruct);
p.addRequired('stimulus',@isstruct);
p.addParameter('bFitZero',@isboolean);
p.parse(params,stimulus,varargin{:});

if p.Results.bFitZero == 0
    % Remove stimulus values of 0 (some analyses might want to leave them in)
    stimulus.stimValues = stimulus.stimValues(stimulus.stimValues~=0);
end

% number of types of parameters, e.g. amp, tau2
numParamTypes = length(params.paramNameCell);

% get unique stimulus values
[uniqueStimValues,~] = unique(stimulus.stimValues);

% for each parameter type and stimulus value
for i = 1:numParamTypes
    for j = 1:length(uniqueStimValues)
       % go into the parameter matrix, find the appropriate column, and
       % pull out all positions with a given unique stim value
       meanParamValues(i,j) = mean(params.paramMainMatrix(uniqueStimValues(j)==stimulus.stimValues,i));
    end    
end

h = figure;
set(h,'Position',[398 380 1099 664]);
for i = 1:numParamTypes
    subplot(1,numParamTypes,i)
    plot(uniqueStimValues,meanParamValues(i,:),'-k^','MarkerSize',10,'LineWidth',1.5); axis square;
    set(gca,'Xscale','log'); set(gca,'Xtick',uniqueStimValues);
    ylabel(params.paramNameCell(i)); xlabel('Temporal frequency');
    set(gca,'FontSize',15);
end

end