# TPUPToolbox

This toolbox implements a three-component model for the fitting of pupil responses to
 a step function of stimulation.

The model assumes (or is optimized for):
	- the stimulus onset begins at the initial time point
	- the response vector is set to have a value of zero at the initial time point
	- the response has largely negative evoked values (a constrictive pupil response)

The fit searches over 6 parameters, three temporal and three gain.
