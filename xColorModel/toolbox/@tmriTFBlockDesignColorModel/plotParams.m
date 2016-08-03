function h = plotParams(obj,params,stimulus,varargin)
% plotParams(obj,param,stimulus,varargin)
%
% Generates plots of parameters for TF Block Design Color Model
%
% Inputs
%   params   : parameter structure
%   stimulus : stimulus structure

p = inputParser;
p.addRequired('param',@isstruct);
p.addRequired('stimulus',@isstruct);
p.parse(params,stimulus,varargin{:});

% Remove stimulus values of 0 (some analyses might want to leave them in)
stimulus.stimValues = stimulus.stimValues(stimulus.stimValues~=0);

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