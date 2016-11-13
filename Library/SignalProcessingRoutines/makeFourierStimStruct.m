function [ stimulusStructOut, fourierSetStructure ] = makeFourierStimStruct( stimulusTimebase, eventTimesArray, msecsToModel, numFourierComponents )
% function [ stimulusStructOut, fourierSetStructure ] = makeFourierStimStruct( stimulusTimebase, eventTimesArray, msecsToModel, numFourierComponents )
%
% This function returns two components that are needed for performing a
% Fourier basis set analysis of events. The resulting stimulusStructureOut
% may be used as the stimulusStructure to model events using the IAMP model
% (with linear regression used to efficiently derive the amplitudes of
% the Fourier components). The fourierSetStructure is combined with the
% amplitudes estimted by the IAMP to produce an estimate of the response.
%
% Arguments in:
%   - stimulusTimebase: the timebase upon which the fourier model will be
%     built (units in msecs)
%   - eventTimesArray: an array of integer values giving the times (in
%     msecs) of events that are then to be modeled with the Fourier basis.
%   - msecsToModel: the duration (in msecs) of the window over which the
%     Fourier basis set will model events. A typical BOLD fMRI application
%     to derive the shape of the HRF would use a window of 16000 msecs.
%   - numFourierComponents: How many sines and cosines (including the
%     zeroeth or dc frequency) to be included in the model. Typically, the
%     number of components to include is equal to the msecsToModel divided
%     by the resolution to be modeled (e.g., use 16 FourierComponents to
%     model the HRF with 1 second resolution over a 16000 msec window).
%
% Arguments out:
%   - stimulusStructOut: a stimulusStruct that has the timebase passed in
%     as stimulusTimeBase, and an m x n matrix of stimulus values, where m
%     is the number of Fourier components and n is the length of the
%     timebase.
%   - fourierSetStructure: a stimulusStruct that has a timebase with the
%     same deltaT as that of stimulusTimebase, and a duration of
%     msecsToModel, and an m x p matrix of stimulus values, where m is the
%     number of Fourier components, and p is the length of the timebase. 
%   

% derive the deltaT from the stimulusTimebase
check = diff(stimulusTimebase);
deltaT = check(1);

% create the fourier set matrix
fourierSetStructure.timebase = ...
    linspace(0,msecsToModel-deltaT,msecsToModel);
componentIndex = 1;
fourierSetStructure.values(1,:) = ...
    fourierSetStructure.timebase*0+1; % Create DC component
componentIndex = componentIndex+1;

% loop through the requested harmonics (numFourierComponents / 2)
for i = 1:ceil(numFourierComponents/2)
    % Create sine for the harmonic
    fourierSetStructure.values(componentIndex,:) = ...
        sin(fourierSetStructure.timebase/msecsToModel*2*pi*i);
    componentIndex = componentIndex+1;
    % Create cosine for the harmonic
    fourierSetStructure.values(componentIndex,:) = ...
        cos(fourierSetStructure.timebase/msecsToModel*2*pi*i);
    componentIndex = componentIndex+1;
end

% trim the set down to the requested number of frequencies
fourierSetStructure.values=fourierSetStructure.values(1:numFourierComponents,:);

% the stimulus values are built first as the vector of attention events
stimulusStructOut.timebase=stimulusTimebase;
impulseEvents=stimulusTimebase*0;
impulseEvents(eventTimesArray)=1;
stimulusStructOut.values = ...
    repmat(impulseEvents,numFourierComponents,1);

% convolve the rows of the attentionEventPacket.stimulus.values by each of
% the rows of the fourierSet
for ii=1:numFourierComponents
    % convolve
    valuesRowConv = conv(stimulusStructOut.values(ii,:), ...
        fourierSetStructure.values(ii,:),'full');
    % cut off extra conv values
    stimulusStructOut.values(ii,:) = ...
        valuesRowConv(1:length(stimulusTimebase));  
end % loop through rows of the Fourier Set

% mean center the component that models the zeroeth frequency events.
stimulusStructOut.values(1,:) = ...
    stimulusStructOut.values(1,:) - ...
    mean(stimulusStructOut.values(1,:));

end

