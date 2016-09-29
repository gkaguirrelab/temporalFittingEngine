function paramsSorted = paramSort4crossVal(obj,params,uniqueStimValues,eventValues)

% function paramsSorted = paramSort4crossVal(obj,paramsVec,uniqueStimValues,eventValues)
%
% NOTE: IF THIS FUNCTION BREAKS, THERE IS A GOOD CHANCE TRANSPOSING ONE OF
% THE INPUTS WILL FIX IT
%
% once you have fit a set of parameters to BDCM data and want to
% cross-validate, this function will be useful. The best way to illustrate
% this is with an example. Say you have three stimulus types, A, B, and C,
% and you have fit a parameter value of 57 to A, a parameter value of 66 to
% B, and a parameter value of 35 to C. In the cross-validation, you are
% confronted with the following sequence of events: BBACABC. This function
% will assign the corresponding sequence of parameters, i.e. [66 66 57 35 57 66 35]
%
% inputs
%
% params          : matrix of parameters. Number of rows = number of unique stimulus
%                   values. Number of columns = number of parameter types.
% uniqueStimValues: unique stimulus values. Must be 1xn vector, where n is
%                   the number of unique stimuli
% eventValues      : stimulus values corresponding to actual events
%
% output
%
% paramsSorted    : vector in which each event has been matched to the
%                   corresponding parameter

p = inputParser;
p.addRequired('params',@isnumeric);
p.addRequired('uniqueStimValues',@isnumeric);
p.addRequired('eventValues',@isnumeric);
p.parse(params,uniqueStimValues,eventValues);

if size(uniqueStimValues,1) ~= 1
   error('paramSort4crossVal: uniqueStimValues MUST be a 1xn vector'); 
end

% create a table matching each stimulus value to the fit parameter
paramStimValueLookup = [uniqueStimValues' params'];

[~,LocInd] = ismember(eventValues,paramStimValueLookup(:,1));

paramsSorted = paramStimValueLookup(LocInd,2:size(paramStimValueLookup,2));

gribble = 1;