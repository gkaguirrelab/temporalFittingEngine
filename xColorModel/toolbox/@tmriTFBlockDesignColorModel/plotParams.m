function h = plotParams(obj,param,stimValues)

% function plotParams(obj,param,stimValues)
%
% generates plots of parameters for TF Block Design Color Model
%
% param     : parameter struct
% stimValues: for a given run, what was the order of stimulus values?

p = inputParser;
p.addRequired('param',@isstruct);
p.addRequired('stimValues',@isnumeric);
p.parse(param,stimValues);

% grab remove stimulus values of 0 (some analyses might want to leave them
% in)
stimValues = stimValues(stimValues~=0);

% number of types of parameters, e.g. amp, tau2
numParamTypes = length(param.paramNameCell);

% get unique stimulus values
[uniqueStimValues,~] = unique(stimValues);

% for each parameter type and stimulus value
for i = 1:numParamTypes
    for j = 1:length(uniqueStimValues)
       % go into the parameter matrix, find the appropriate column, and
       % pull out all positions with a given unique stim value
       meanParamValues(i,j) = mean(param.paramMainMatrix(uniqueStimValues(j)==stimValues,i));
    end    
end

h = figure;
set(h,'Position',[398 380 1099 664]);
for i = 1:numParamTypes
    subplot(1,numParamTypes,i)
    plot(uniqueStimValues,meanParamValues(i,:),'-k^','MarkerSize',10,'LineWidth',1.5); axis square;
    set(gca,'Xscale','log'); set(gca,'Xtick',uniqueStimValues);
    ylabel(param.paramNameCell(i)); xlabel('Temporal frequency');
    set(gca,'FontSize',15);
end

end