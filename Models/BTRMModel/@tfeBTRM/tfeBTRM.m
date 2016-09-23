classdef tfeBTRM < tfe
% tfeBTRM
%
%   tfe =  tfeBTRM(vargin);
% 
% Implements a model for a block design fMRI experiment in which the color
% direction and temporal frequency of stimuli are varied in a block design.
%
% Inherits optional key/value pairs from parent class tfe.
%
% 7/25/16  dhb, bc  Started on this.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        paramsBase;
    end
    
    % Public methods
    methods  
    end
    
    properties (Dependent)
    end
    
    % Methods.  Most public methods are implemented in a separate function,
    % but we put the constructor here.
    methods (Access=public)
        % Constructor
        function obj = tfeBTRM(varargin)
            
            % Base class constructor
            obj = obj@tfe(varargin{:});
            
            % You could do model specific stuff here.

        end
    end 
    
    % Get methods for dependent properties
    methods
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
