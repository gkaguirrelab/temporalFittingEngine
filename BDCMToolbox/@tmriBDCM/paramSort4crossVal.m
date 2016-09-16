function paramsSorted = paramSort4crossVal(obj,params,uniqueStimValues,stimValues)

% function paramsSorted = paramSort4crossVal(obj,paramsVec,uniqueStimValues,stimValues)
%
% params          : matrix of parameters. Number of rows = number of unique stimulus
%                   values. Number of columns = number of parameter types.
% uniqueStimValues: unique stimulus values. Must be 1xn vector, where n is
%                   the number of unique stimuli
% stimValues      : stimulus values corresponding to events

p = inputParser;
p.addRequired('paramsVec',@isnumeric);
p.addRequired('uniqueStimValues',@isnumeric);
p.addRequired('stimValues',@isnumeric);
p.parse(params,uniqueStimValues,stimValues);

if size(uniqueStimValues,1) ~= 1
   error('paramSort4crossVal: uniqueStimValues MUST be a 1xn vector'); 
end

paramStimValueLookup = [params' uniqueStimValues'];

gribble = 1;