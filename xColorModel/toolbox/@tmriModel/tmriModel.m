classdef tmriModel < handle
% The @tmri parent class for modeling temporal responses.
% 
% 6/26/16  dhb  Started in on this

    % Public read/write properties
    properties
        % Structure of parameters for the current simulated model. The
        % user can set these and then get the simulated model response.
        simulateParams = [];
        
        % Timebase on which to compute model predictions, specified in
        % seconds.
        timebase = [];
        
        % Stimulus.  The exact form this takes is model dependent, but it
        % should allow the compute methods to predict the stimulus on the
        % timebase, given the parameters of the model.
        stimulus = [];
            
        % HRF
        boldHRF = [];
  
    end
    
    % Dependent properties, computed from other parameters
    properties
        % Neural response on timebase
        simulateNeuralResponse = [];
        fitNeuralResponse = [];
        
        % Bold response
        simulateBoldResonse = [];
        fitBoldResponse = [];
    end
        
    % Public, read-only properties.  These can be set by methods of the
    % parent class (that is, this class) but not by methods of subclasses.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected, GetAccess = public)
        % Structure of parameters for the current fit.
        fitParams = [];

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
        % Default parameters
        defaultParams(obj,varargin);
        
        % Print parameters
        printParams(obj,varagin);

        % Convert parameter struct to a vector to be used by search
        % routines.
        x = paramsToVec(obj);
        
        % Take the vector and put it back into the object's parameter
        % structure.
        vecToParams(obj,x)
        
        % Compute forward simulation of the implemented model, given the parameters
        computeNeural(obj,varargin);
          
        % Might want the object to know how to plot itself.
        %tmriPlot(obj, sensor);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
