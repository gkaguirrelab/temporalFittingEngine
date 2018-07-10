TPUP (three-component pupil) model.

Models an evoked pupil response with a 6-parameter, 3-component model.

The input to the model is the stimulus profile. An additional two input
 vectors, representing the rate of stimulus change at onset, are created
 by differentiating the stimulus profile and retaining the positive
 elements. These three vectors are then subjected to convolution
 operations composed of a gamma and exponential decay function, each
 under the control of a single time-constant parameter. The resulting
 three components (red) were normalized to have unit area, and then
 subjected to multiplicative scaling by a gain parameter applied to each
 component. The scaled components are summed to produce the modeled
 response, which is temporally shifted.

The response to be modeled should be in % change units (e.g. 10%
 contraction, as opposed to 0.1) so that the various parameters have
 similar magnitudes of effect upon the modeled response.

delay - time to shift the model to the right (msecs)

gammaTau - time constant of the Gamma function (msecs)

exponentialTau - time constant of the persistent component (seconds)

amplitudeTransient - scaling of the transient component in (%change*secs)

amplitudeSustained - scaling of the transient component in (%change*secs)

amplitudePersistent - scaling of the transient component in (%change*secs)
