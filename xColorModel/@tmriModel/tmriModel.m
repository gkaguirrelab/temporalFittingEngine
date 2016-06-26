classdef tmriModel < handle
% The @tmri parent class for modeling temporal responses.
% 
% 6/26/16  dhb  Started in on this
% 
    % Public read/write properties
    properties
        % Structure of parameters for the current model
        params;
              
        % Timebase on which to compute model predictions
        timebase;
        
        % Stimulus.  The exact form this takes is model dependent,
        % but it should be a vector specified on the timebase.
        stimulus;

        % Noise flag.  Add noise to forward simulated data?
        noiseflag;
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        % Prediction. The model's prediction on the timebase, given
        % the current parameters.
        prediction;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    % Public methods
    methods
        function obj = tmriModel(varargin)
            % obj.initialize(varargin{:});
            obj.params = [];
            obj.timebase = [];
            obj.stimulus = [];
            obj.noiseflag = false;

        end
        
        % Call set through a setter function.  I don't think we need this.
        % function obj = set(obj, varargin)
        %     osSet(obj, varargin{:});
        % end
       
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % Compute forward simulation of the implemented model, given the parameters
        compute(obj,varargin);
        
        % Convert parameter struct to a vector to be used by search
        % routines.
        x = paramstovec(obj);
        
        % Take the vector and put it back into the object's parameter
        % structure.
        vectoparams(obj,x)
        
        % Might want the object to know how to plot itself.
        %tmriPlot(obj, sensor);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        %initialize(obj);
    end
    
end
