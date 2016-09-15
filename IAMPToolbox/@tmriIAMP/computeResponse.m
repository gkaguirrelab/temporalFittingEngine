function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the instance amplitude model.
%
% The respone is simply the sum of the stimulus vectors, after each
% stimulus vector has been multiplied by an amplitude parameter for that
% instance.
%
% Optional key/value pairs
%   'AddNoise'
%     true or false(default) 
%  'HRF' - a structure describing a kernel to be used to go from neural to observed response.
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
% *assume timebase is the same for all stimuli*
individualResponses = forwardModelIAMP(timebase, stimulus.values,...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'Amplitude')));
summedResponse = sum(individualResponses,1);

%% Optionally, convolve with a kernel
response = obj.applyHRF(timebase,summedResponse,p.Results.HRF);

%% Optional add noise
if (p.Results.AddNoise)
    response = response + normrnd(0,params.noiseSd,size(response));
end

end