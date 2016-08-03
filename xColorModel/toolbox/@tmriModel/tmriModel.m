classdef tmriModel < handle
% The @tmri parent class for modeling temporal responses.
% 
% 6/26/16  dhb  Started in on this

    % Public read/write properties
    properties              
    end
    
    % Dependent properties, computed from other parameters
    properties
    end
        
    % Public, read-only properties.  These can be set by methods of the
    % parent class (that is, this class) but not by methods of subclasses.
    properties (SetAccess = private, GetAccess = public)  
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set or read these
    properties (Access = private)      
    end
    
    % Public methods
    %
    % Methods defined in separate files are public by default, so we don't
    % explicitly decare them.  But, if you wanted to write the whole body
    % of some short public method here, you could do so.
    methods (Access=public)
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % Default parameters and bounds in structure form
        [params,paramsLb,paramsUb] = defaultParams(obj,varargin);
        
        % Print parameters
        print(obj,params,varagin);
        
        % Get parameter locking matrix for this model
        paramLockMatrix = lockMatrix(obj,params,varargin);

        % Convert parameter struct to a vector to be used by search
        % routines.
        x = paramsToVec(obj,params);
        
        % Take the vector and put it back into the object's parameter
        % structure.
        params = vecToParams(obj,x)
        
        % Compute forward simulation of the implemented model, given the parameters
        response = computeResponse(obj,params,timebase,stimulus,varargin);
        
        % Plot fit parameters for neural model
        h = plotParams(obj,param,stimValues);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
