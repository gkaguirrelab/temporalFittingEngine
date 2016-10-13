# FCON -- Effective Contrast Model

This model acts as a wrapper to expand a single parameter (effective contrast) into a set of parameters, which then call the forward model of a passed tfe model subclass. This allows prior knowledge of the possible patterns in parameters to constrain the available fits.

A key aspect of the use of the FCON model is the generation of an fcon struct.  The fcon struct is passed within packet.stimulus, and has the fields:

    fcon.contrastbase - the n contrast levels at which an effective contrast value is available
    fcon.paramLookUpMatrix - a m x n matrix where m is the number of parameters in the expanded description of the data, and n is the
       number of effective contrast levels.
    fcon.modelObjHandle - an object handle that identifies the model subclass to be used to calculate the forward model using the expanded parameter set.
    fcon.defaultParams - the default params for the model subclass

The contrastbase and paramLookUpMatrix is a look-up table for relating an effective contrast value to an expanded set of parameters.

The fmincon search is performed by gradient descent. It will become stuck if small changes in the parameter (effective contrast) do not produce any change in the objective function. Because we are searching not over a continuous function but instead through a discontinuous look-up table of mappings between effective contrast and expanded parameters, we need to know what is the size of the spacing between adjacent effective contrast values, and tell fmincon to have the minimum gradient step size in the search be at least this big. Therefore, it is recommended that:

  - The spacing between values in the contrastbase be regular. For some applications, it may be desirable to have the contrastbase vector actually contain log10 of the contrast values.
  - The key-value pair 'DiffMinChange' should be set equal to the spacing between values in contrastbase, and passed to the tfe when calling the fitResponse method
  
 