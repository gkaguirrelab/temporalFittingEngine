function h = plotParams(obj,params,stimulus,varargin)
% plotParams(obj,param,stimulus,varargin)
%
% Generates plots of parameters.  Currently a stub.
%
% Inputs
%   params   : parameter structure
%   stimulus : stimulus structure

p = inputParser;
p.addRequired('param',@isstruct);
p.addRequired('stimulus',@isstruct);
p.parse(params,stimulus,varargin{:});

% Here is where one would implement a plot and set the figure handle to it.
h = [];

end