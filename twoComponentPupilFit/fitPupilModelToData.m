function [yPupil, fitParams, fitFigureHandle] = fitPupilModelToData(t, data, stdErrorVector, stimulus, initialParamsIn)

%% fitPupilModelToData
%
% This function fits a two stage, two component pupil response model to
%  pupil response time-series. The model includes:
%
%   - Convolution of the stimulus profile with a gamma function
%   - Sustained response, modeled as a decaying exponential
%   - Persistent response, modeled as a supra-saturating function
%
% Input properties:
%
%   t - a vector of time points.
%   data - a vector of pupil amplitudes
%   stimulus - a vector of stimulus amplitudes
%   initialParamsIn - a vector of parameters (optional)
%
% 07-01-16 -  gka wrote it
% 07-02-16 - refinements to the model (gka)
%

% Assume that we do not want to plot the fit unless we hear otherwise
displayFitPlot = true;
fitFigureHandle = nan;

%% Default parameter definition
param.startTime = 0.1; % time (in seconds) that the initial transient pupil response begins

% parameters of the initial gamma convolution
param.gammaTau = 0.2;    % time constant of the Gamma function

% parameters of the sustained response
param.sustainedAmp = 0.2; % scaling of the sustained component
param.sustainedTau = 0.5;  % time constant of the low-pass (exponential decay) component.

% parameters of the persistent
param.persistentAmp = 0.2; % amplitude scaling of the persistent response
param.persistentT50 = 0.25; % time to half-peak of the super-saturating function
param.persistentAlpha = 1.5;  % time constant of the decay of the super-saturating function.

% convert the error vector to serve as a weight vector
weightVector=1./stdErrorVector.^2;

if nargin==5
    initialParams=initialParamsIn;
else
    initialParams(1) = param.startTime;
    initialParams(2) = param.gammaTau;
    initialParams(3) = param.sustainedAmp;
    initialParams(4) = param.sustainedTau;
    initialParams(5) = param.persistentAmp;
    initialParams(6) = param.persistentT50;
    initialParams(7) = param.persistentAlpha;
end

% Set some bounds
LowerBounds=[0.1,0.1,0,0,0,0,0.5];
UpperBounds=[0.4,0.3,1,2,1,2,6.0];

% Some of the pupil data has non-zero values for the response as
% compared to the stationary background. We offset the data to have the
% lowest value at zero, and retain this to adjust the fit later.
offset=nanmean(data(1:100));
data=data-offset;

% Create the optimization options structure, and increase the limit of
% functions to evaluate
options = optimset('MaxFunEvals',1000);
options = optimset('MaxIter',1000);

% Run the optimization
[fitParams,fval,exitflag] = fminsearchbnd( (@(p) pupilFitError(p, t, data, weightVector, stimulus)), initialParams, LowerBounds, UpperBounds, options);

% Obtain the pupil model fit
[yPupil,gammaStimulus,ySustained,yPersistent] = createPupilTemporalModel(t,stimulus,fitParams);

% add in the offset
yPupil=yPupil+offset;
gammaStimulus=gammaStimulus+offset;
ySustained=ySustained+offset;
yPersistent=yPersistent+offset;

% If the user requested a plot, give it to them

if displayFitPlot
    stimulusPlotScale=0.25; % Amount to scale stimulus vectors in plot for display
    fitFigureHandle=figure();
    % Plot the data (with ± error bounds)
    plot(t, data+offset+stdErrorVector, '.', 'color', [0.8 0.8 0.8]);
    hold on;
    r0 = plot(t, data+offset-stdErrorVector, '.', 'color', [0.8 0.8 0.8]);
    r1 = plot(t, data+offset, '.r'); hold on;
    
    % Plot the  fits
    r2 = plot(t, yPupil, '-k', 'LineWidth', 2);
    r3 = plot(t, stimulus*stimulusPlotScale);  % Scale the stimulus plots to keep them on the same axis
    r4 = plot(t, gammaStimulus*stimulusPlotScale);  % Scale the stimulus plots to keep them on the same axis
    r5 = plot(t, ySustained*(-1));
    r6 = plot(t, yPersistent*(-1));
    
    % Make the plot pretty
    xlabel('time');
    ylabel('Amplitude');
    legend([r0 r1 r2 r3 r4 r5 r6], '±Error', 'Data', 'Fit', 'stimulus (scaled/4)', 'gamma Stimulus (scaled/4)', 'Sustained', 'persistent'); legend boxoff;
    hold off;
end

gribble=1;

function E = pupilFitError(params, t, data, weightVector, stimulus)

% Error function, calculating the sum-of-squares for the data vs. the fit.

[model,~,~,~] = createPupilTemporalModel(t, stimulus, params);

% Calculate the sums-of-squares

E = sum(((data - model).^2).*weightVector);


