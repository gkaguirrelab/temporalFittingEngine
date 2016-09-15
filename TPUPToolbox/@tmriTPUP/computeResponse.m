function response = computeResponse(obj,params,timebase,stimulus,varargin)
% response = computeResponse(obj,params,timebase,stimulus,varargin)
% 
% Compute method for the two-component pupil response model.
%
% The forward model consists of two components (an exponentially decaying
%   sustained response, and a super-saturating function persistent
%   response), under the control of 7 parameters.
%
% Optional key/value pairs
%   'AddNoise'
%     true or false(default) 
%  'HRF' - a structure describing a kernel to be used to be applied after
%    the forward model. Empty matrix is default, in which case no
%    convolution is performed done

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


%% Compute the forward model
% *assume timebase is the same for all stimuli*
individualResponses = forwardModel(timebase, stimulus.values,...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'startTime')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'gammaTau')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedAmp')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedTau')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAmp')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'Amplitude')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentT50')), ...
    params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAlpha'))     );
summedResponse = sum(individualResponses,1);

%% Optionally, convolve with a passed kernel
response = obj.applyHRF(timebase,summedResponse,p.Results.HRF);

%% Optional add noise
if (p.Results.AddNoise)
    response = response + normrnd(0,params.noiseSd,size(response));
end

end