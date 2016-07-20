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
    end
    
    properties (Dependent)
        % The quadratic matrix.  This depends on the parameters, but we
        % often want to use it explicity.
        simulateQ
        fitQ
    end
    
    % Methods that must be implemented (Abstract in parent class).
    methods (Access=public)
        % Return reasonable default parameters for the model
        function [paramsVec,vlbVec,vubVec] = defaultParams(obj,varargin)
            [paramsVec,vlbVec,vubVec] = tmriDefaultParams(obj,varargin);
        end
        
        % Print the parameters
        function printParams(obj,varargin)
            tmriPrintParams(obj,varargin);
        end
        
        % Convert parameter struct to a vector to be used by search
        % routines.
        function x = paramsToVec(obj,varargin)
            x = tmriParamsToVec(obj,varargin);
        end
        
        % Take the vector and put it back into the object's parameter
        % structure.
        function vecToParams(obj,x,varargin)
            tmriVecToParams(obj,x),varargin;
        end
                  
        % Forward simulation of the implemented model, given the parameters
        function obj = computeNeural(obj,varargin)
            obj = tmriComputeNeural(obj,varargin);
        end
    end 
    
    % Get methods for dependent properties
    methods
        % Get quadratic form in matrix form from the parameters
        function q = get.simulateQ(obj)
            [~,~,q] = EllipsoidMatricesGenerate([1 obj.simulateParams.Qvec]');
        end
        
        function q = get.fitQ(obj)
            [~,~,q] = EllipsoidMatricesGenerate([1 obj.fitParams.Qvec]');
        end
        
       function response = get.simulateNeuralResponse(obj)
            response = tmriComputeNeural(obj,varargin);
       end
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
