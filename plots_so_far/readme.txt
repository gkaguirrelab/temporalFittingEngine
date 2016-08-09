asb1 and gka1 are subjects

Filenames and meaning:

(subject)_carryover_(modulation direction)_amp: 

carry over matrix for amplitude. Middle panel shows the direct effect of each stimulus (marginalize carry over matrix across columns), and the right panel shows the 'preceding effect' (marginalize over rows)

(subject)_carryover_(modulation direction)_tau2: 

carry over matrix for the exponential decay parameter. Doesn't look like there is much to see here.

(subject)_fits: 

Averaged time series, with fits. For when parameters for all stimuli with the same temporal frequency are locked together. Fits are OK, but not great. Blue line is the data with standard error shading, red line is the fits. Numbers on top indicate the temporal frequency of the stimulus starting at that point in time. Corresponding grayscale bars below indicate the temporal frequency of the stimulus during that interval: white corresponds to 0Hz (no stimulus), black to 64Hz. Mean squared error also indicated. 

(subject)_unlocked: 

fits, with each 12-second block allowed to have its own set of parameters. Good fits. 

(subject)_params:

Plots for the different parameters. 

(subject)_fits_closeup:

Used to examine fit of the model at a high level of detail. For each combination of modulation direction and temporal frequency, find all stimulus onsets. Then define a 24 second window for each onset, and average them. 