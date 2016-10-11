# FCON -- Effective Contrast Model

This model acts as a wrapper to expand a single parameter (effective contrast) into a set of parameters, which then call the forward model of a passed tfe model subclass. This allows prior knowledge of the possible patterns in parameters to constrain the available fits.