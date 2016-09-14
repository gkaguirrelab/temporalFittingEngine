function [response]=createDoubleGammaFunction(t,paramsIn)

% vec to params
    
params_gamma1shape = paramsIn(1); 
params_gamma2shape = paramsIn(2); 
params_gammaScale = paramsIn(3);  % relative amplitude of neg to pos gamma

% Generate the model. We return response.

response =  1 * ... % Amplitude not important
            (gampdf(t, params_gamma1shape, 1) - ...
             gampdf(t, params_gamma2shape, 1)*params_gammaScale);
 
end