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
            
            % Set default parameters
            tmriDefaultParams(obj,varargin);
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
    
    properties (Dependent)
        % The quadratic matrix.  This depends on the parameters, but we
        % often want to use it explicity.
        Q
    end
    
    % Methods that must be implemented (Abstract in parent class).
    methods (Access=public)
        % Return reasonable default parameters for the model
        function [paramsVec,vlbVec,vubVec] = defaultParams(obj)
            [paramsVec,vlbVec,vubVec] = tmriDefaultParams(obj);
        end
        
        % Convert parameter struct to a vector to be used by search
        % routines.
        function x = paramsToVec(obj)
            x = tmriParamsToVec(obj);
        end
        
        % Take the vector and put it back into the object's parameter
        % structure.
        function vecToParams(obj,x)
            tmriVecToParams(obj,x);
        end
                  
        % Forward simulation of the implemented model, given the parameters
        function obj = computeNeural(obj,varargin)
            obj = tmriComputeNeural(obj,varargin);
        end
    end 
    
    % Get methods for dependent properties
    methods
        % Get quadratic form in matrix form from the parameters
        function q = get.Q(obj)
            [~,~,q] = EllipsoidMatricesGenerate([1 ; obj.params.Qvec]);
        end
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
