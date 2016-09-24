function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the two-component pupil response model.
%
% Operates by calling forwardModelTPUP and then the applyKernel tfe method.
%
% Optional key/value pairs
%   'AddNoise' - true/false (default false).  Add noise to computed
%     response?  Useful for simulations.

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('stimulusStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) || isstruct(x)));
p.addParameter('addNoise',false,@islogical);
p.parse(params,stimulusStruct,kernelStruct,varargin{:});
params = p.Results.params;

%% Detect the presence of NaN values in the paramMainMatrix.
% Not sure why this happens sometimes. If we find NaNs, return
% a response vector of zeros, and update the screen with an 'n'
if sum(isnan(params.paramMainMatrix))~=0
    fprintf('n');
    response=timebase*0;
else
    
    %% Compute the forward model
    % *assume timebase is the same for all stimuli*
    % GEOFF PICK IT UP HERE.  NOTE FROM DHB: I UPDATED THE CONVOLUTION AND
    % NOISE ADDITION TO MATCH WHAT I THOUGHT IT WOULD EVENTUALLY BE, BASED
    % ON WHAT I WAS DOING IN THE QCM VERSION.
    individualResponses = forwardModelTPUP(params,stimulusStruct,...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'startTime')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'gammaTau')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedAmp')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'sustainedTau')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAmp')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentT50')), ...
        params.paramMainMatrix(:,strcmp(params.paramNameCell,'persistentAlpha'))     );
    response = sum(individualResponses,1);
    
    % report an iteration has completed
    fprintf('.');
    
    %% Optionally, convolve with a passed kernel
    modelResponseStruct = obj.applyKernel(modelResponsStruct,kernelStruct,varargin{:});
    
    %% Optional add noise
    if (p.Results.addNoise)
        modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(response));
    end
end

end