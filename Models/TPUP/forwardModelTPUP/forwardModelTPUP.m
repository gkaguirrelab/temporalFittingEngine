function [modelResponseStruct] = forwardModelTPUP(obj,params,stimulusStruct)

%% forwardModelTPUP
%
% Models the pupil temporal response as three, temporally
% overlapping components, each controlled with an amplitude and one or
% two time-constant parameters.
%
%  The passed stimulus vector is first differentiated to find the onset.
%    The convolved with a gamma function,
%    with one parameter that defines the time constant (shape). This
%    provides some of the "peakiness" of the initial transient seen in the
%    response data. There are then two components:

%  Transient -- A gamma function controlled by a shape and gain parameter,
%     positioned at the time of onset of the events
%  Sustain -- The stimulus subjected to an exponential decay
%  Persistent -- The convolved stimulus vector is subjected to
%     convolution with a super-saturating function
%
% An additional parameter allows the entire model to shift forward or back
% in time relative to the data.
%



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
    
    transientComponent = obj.applyKernel(stimulusSlewOn,gammaIRF);
    sustainedComponent = obj.applyKernel(stimulus,gammaIRF);
    persistentComponent = obj.applyKernel(obj.applyKernel(stimulusSlewOn,exponentialIRF),gammaIRF);
    
    % Scale each component to have unit area
    transientComponent=normalizeKernelArea(transientComponent);
    sustainedComponent=normalizeKernelArea(sustainedComponent);
    persistentComponent=normalizeKernelArea(persistentComponent);
    
    yPupil=transientComponent.values * amplitudeTransietVec(ii) + ...
        sustainedComponent.values * amplitudeSustainedVec(ii) + ...
        persistentComponent.values * amplitudePersistentVec(ii);
    
    % apply the temporal delay
    yPupil=fshift(yPupil,-1*delayVec(ii)/deltaT);
    yPupil(1:ceil(delayVec(ii)/deltaT))=0;
    
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
