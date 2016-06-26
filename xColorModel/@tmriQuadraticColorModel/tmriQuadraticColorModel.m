classdef tmriQuadraticColorModel < tmriModel
% tmriQuadraticColorModel
%
%   tmri = tmriQuadraticColorModel();
% 
% Implements a model that is quadratic in the color contrast of the
% stimulus.
%
% 6/26/16  dhb  Started in on this.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % You can define some properties here.
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = tmriQuadraticColorModel(varargin)
            % Initialize the parent class
            obj = obj@tmriModel();
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must be implemented (Abstract in parent class).
    methods (Access=public)
        
        % Forward simulation of the implemented model, given the parameters
        function obj = compute(obj,varargin)
            obj = tmriCompute(obj,varargin);
        end
        
        % Convert parameter struct to a vector to be used by search
        % routines.
        function x = paramstovec(obj)
            x = tmriParamstovec(obj);
        end
        
        % Take the vector and put it back into the object's parameter
        % structure.
        function vectoparams(obj,x)
            obj = tmriVectoparams(obj,x);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
