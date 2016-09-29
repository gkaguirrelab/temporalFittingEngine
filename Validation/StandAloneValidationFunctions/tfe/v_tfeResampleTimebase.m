% resampleTimebaseValidate
%
% Assessments of the resampleTimebase method in the parent tfe class

%% Clear and close
clear; close all;

%% Construct the model object
% We have to pick some subclass to define the model, although this
% shouldn't matter for the validation tests here
temporalFit = tfeIAMP('verbosity','high');

%% Perform some tests

% Test if brief events at high temporal sampling are represented after
% down-sampling
deltaT=10; % msecs
totalTime=3000; % msecs
impulseLocation=1500; % msecs

inputStruct.timebase=linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(inputStruct.timebase,2);
inputStruct.values=zeros(1,nTimeSamples);
inputStruct.values(1,round(impulseLocation/deltaT))=1;

temporalFit.plot(inputStruct,'Marker','x')

% create the new timebase
newDeltaT=100; % msecs
newTimebase=linspace(1,totalTime,newDeltaT);

% obtain the resampled testStruct
outputStruct=temporalFit.resampleTimebase(inputStruct,newTimebase);

temporalFit.plot(outputStruct,'Marker','x','Color',[0 1 0])

% Assess the resampled plot
expectedMax=max(inputStruct.values)*(deltaT/newDeltaT);
