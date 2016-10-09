# IAMPToolbox

This is an elementary model implemented within the temporal fitting engine.
It simply estimates an amplitude for each stimulus instance (then convolved by a kernel)
to fit the data vector.

The linearRegression search method is available for this model, which is much faster than non-linear (fmincon) searches.