function [response]=createSuperSaturatingFunction(t,paramsIn)

% Credit to n.c. for this function.

    % vec to params
    
    param.t_50=paramsIn(1); % time (in seconds) of the peak
    param.alpha=paramsIn(2); % time constant of the fall post peak

    % fixed parameters of the function
    
    n=2;
    Ro=0;
    Rmax=1.0;
    
    response = Rmax * (t.^n)./(t.^(param.alpha*n) + param.t_50^(param.alpha*n)) + Ro;
end