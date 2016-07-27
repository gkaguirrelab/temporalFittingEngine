function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the BDCM model.
%
% Optional key/value pairs
%   'AddNoise'
%     true or false(default) 
%  'HRF' - a structure describing the HRF to be used to go from neural to BOLD response.
%    Empty matrix is default, in which case no convolution is done

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('timebase',@isnumeric);
p.addRequired('stimulus',@isstruct);
p.addParameter('AddNoise',false,@islogical);
p.addParameter('HRF',[]);
p.parse(params,timebase,stimulus,varargin{:});
params = p.Results.params;
timebase = p.Results.timebase;
stimulus = p.Results.stimulus;

%% First compute the neural response
individualResponses = createNeuralTemporalModelFromStimMatrix(timebase, stimulus.stimMatrix,...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'Amplitude')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'tau2')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'ARAmplitude')));
neuralResponse = sum(individualResponses,1);

%% Optionally, apply HRF
response = obj.applyHRF(timebase,neuralResponse,p.Results.HRF);

%% Optionally, apply HRF
response = obj.applyHRF(timebase,response,p.Results.HRF);

%% Optional add of noise
if (p.Results.AddNoise)
    response = response + normrnd(0,params.noiseSd,size(response));
end

end