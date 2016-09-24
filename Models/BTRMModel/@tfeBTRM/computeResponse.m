function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the BTRM model.
%
% Optional key/value pairs
%   'AddNoise' - true/false (default false).  Add noise to computed
%     response?  Useful for simulations.
%  'HRF' - structure (default empty).  Structure describing HRF to be used
%     to go from neural to BOLD response. If empty, no convolution is done

%% Parse input.
% At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('timebase',@isnumeric);
p.addRequired('stimulus',@isstruct);
p.addParameter('addNoise',false,@islogical);
p.addParameter('HRF',[],@isstruct);
p.parse(params,timebase,stimulus,varargin{:});
params = p.results.params;
timebase = p.results.timebase;
stimulus = p.results.stimulus;

%% First compute the neural response
% *assume timebase is the same for all stimuli*
individualResponses = createNeuralTemporalModelFromStimMatrix(timebase, stimulus.values,...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitude')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'tau2')));
neuralResponse = sum(individualResponses,1);

%% Optionally, apply HRF
response = obj.applyHRF(timebase,neuralResponse,p.results.HRF);

%% Optional add of noise
if (p.results.addNoise)
    response = response + normrnd(0,params.noiseSD,size(response));
end

end