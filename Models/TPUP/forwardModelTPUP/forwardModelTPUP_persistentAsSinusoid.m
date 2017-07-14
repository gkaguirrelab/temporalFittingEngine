function [modelResponseStruct] = forwardModelTPUP(obj,params,stimulusStruct)

%% forwardModelTPUP
%
% Models an evoked pupil response with a 6-parameter, 3-component model.
%
% The input to the model is the stimulus profile. An additional two input
%  vectors, representing the rate of stimulus change at onset, are created
%  by differentiating the stimulus profile and retaining the positive
%  elements. These three vectors are then subjected to convolution
%  operations composed of a gamma and exponential decay function, each
%  under the control of a single time-constant parameter. The resulting
%  three components (red) were normalized to have unit area, and then
%  subjected to multiplicative scaling by a gain parameter applied to each
%  component. The scaled components are summed to produce the modeled
%  response, which is temporally shifted.
%
% The response to be modeled should be in % change units (e.g. 10%
%  contraction, as opposed to 0.1) so that the various parameters have
%  similar magnitudes of effect upon the modeled response.
%
% delay - time to shift the model to the right (msecs)
% gammaTau - time constant of the Gamma function (msecs)
% exponentialTau - time constant of the persistent component (seconds)
% amplitudeTransiet - scaling of the transient component in (%change*secs)
% amplitudeSustained - scaling of the transient component in (%change*secs)
% amplitudePersistent - scaling of the transient component in (%change*secs)



%% Unpack the params
%   Overall model timing
%     startTime -  left or right tie shift for the entire model, relative
%     to the stimulus
%     gammaTau - time constant of the transient gamma function

delayVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'delay'));
gammaTauVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'gammaTau'));
exponentialTauVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'exponentialTau')).*1000;
amplitudeTransietVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitudeTransiet')).*1000;
amplitudeSustainedVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitudeSustained')).*1000;
amplitudePersistentVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'amplitudePersistent')).*1000;
% new sinusoid parameters
% sinusoidCycleTimeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sinusoidCycleTime'));
sinusoidCycleTimeVec=params.paramMainMatrix(:,strcmp(params.paramNameCell,'sinusoidCycleTime')).*1000;

%sinusoidCycleTimeVec=12000; %12,000 ms being a reasonable estimate from looking at the group average

% derive some basic properties of the stimulus values
numInstances=size(stimulusStruct.values,1);
check = diff(stimulusStruct.timebase);
deltaT = check(1);


% pre-allocate the responseMatrix and kernels here for speed
responseMatrix=zeros(numInstances,length(stimulusStruct.timebase));
gammaIRF.values = stimulusStruct.timebase .* 0;
gammaIRF.timebase = stimulusStruct.timebase;
exponentialIRF.values=stimulusStruct.timebase .* 0;
exponentialIRF.timebase=stimulusStruct.timebase;
sinusoidIRF.values=stimulusStruct.timebase .* 0;
sinusoidIRF.timebase=stimulusStruct.timebase;


%% We loop through each row of the stimulus matrix
for ii=1:numInstances
    
    % grab the current stimulus
    stimulus.values = stimulusStruct.values(ii,:);
    stimulus.timebase = stimulusStruct.timebase;
        
    % Create a stimulusSlewOn vector. This is the rate of change of the
    % stimulus at the time of onset.
    stimulusSlewOn.values= max( [ [diff(stimulus.values) 0]; zeros(1,length(stimulus.timebase)) ] );
    stimulusSlewOn.timebase=stimulus.timebase;

    % Create the gamma kernel
    gammaIRF.values = stimulus.timebase .* exp(-stimulus.timebase./gammaTauVec(ii));
    gammaIRF=normalizeKernelArea(gammaIRF);
    
    % Create the exponential kernel
    exponentialIRF.values=exp(-1/exponentialTauVec(ii)*stimulus.timebase);
    exponentialIRF=normalizeKernelArea(exponentialIRF);
    
    % Create sinusoid kernel
    sinusoidIRF.values=sin(stimulus.timebase/sinusoidCycleTimeVec(ii)*2*pi);
    sinusoidIRF.values(round(sinusoidCycleTimeVec(ii)/2/deltaT):length(sinusoidIRF.values)) = 0;
    sinusoidIRF=normalizeKernelArea(sinusoidIRF);
    
    transientComponent = obj.applyKernel(stimulusSlewOn,gammaIRF);
    sustainedComponent = obj.applyKernel(stimulus,gammaIRF);
    %persistentComponent = obj.applyKernel(obj.applyKernel(stimulusSlewOn,exponentialIRF),gammaIRF);
    persistentComponent = obj.applyKernel(stimulus,sinusoidIRF);
    persistentComponent = obj.applyKernel(stimulusSlewOn, stepPulseStruct)
    
    %playing with the sustained component was well
    %sustainedComponent = obj.applyKernel(obj.applyKernel(stimulus,exponentialIRF),gammaIRF);
    
    % Scale each component to have unit area
    transientComponent=normalizeKernelArea(transientComponent);
    sustainedComponent=normalizeKernelArea(sustainedComponent);
    persistentComponent=normalizeKernelArea(persistentComponent);
    
    yPupil=transientComponent.values * amplitudeTransietVec(ii) + ...
        sustainedComponent.values * amplitudeSustainedVec(ii) + ...
        persistentComponent.values * amplitudePersistentVec(ii);
    
    % apply the temporal delay
    initialValue=yPupil(1);
    yPupil=fshift(yPupil,-1*delayVec(ii)/deltaT);
    yPupil(1:ceil(-1*delayVec(ii)/deltaT))=initialValue;
    
    % Add this stimulus model to the response matrix
    responseMatrix(ii,:)=yPupil;
    
end % loop over stimulus instances

% Check the result for nans
if ~sum(sum(isnan(responseMatrix)))==0
    error('NaNs detected in the responseMatrix');
end

%% Build the modelResponseStruct to return
modelResponseStruct.timebase=stimulusStruct.timebase;
modelResponseStruct.values=sum(responseMatrix,1);

end % function
