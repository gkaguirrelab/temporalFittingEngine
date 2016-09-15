# TPUPToolbox

This toolbox implements a two-component model for the fitting of pupil responses to
 a step function of stimulation.

The model assumes (or is optimized for):
	- the stimulus onset begins at the initial time point
	- the response vector is set to have a value of zero at the initial time point
	- the response has largely negative evoked values (a constrictive pupil response)

The fit searches over 7 parameters:
% initialValue - value at the first time point (measured directly from the response values)
% startTime - time (in seconds) that the initial transient pupil response begins
% gammaTau - time constant of the Gamma function (msecs)
% sustainedAmp - scaling of the sustained component
% sustainedTau - time constant of the low-pass (exponential decay) component (seconds) of the sustained response
% persistentAmp - amplitude scaling of the persistent response
% persistentT50 - time to half-peak of the super-saturating function (seconds)
% persistentAlpha - time constant of the decay of the super-saturating function (seconds).
