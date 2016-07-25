classdef tmriTFBlockDesignColorModel < tmriModel
% tmriQuadraticColorModel
%
%   tmri =  tmriTFBlockDesignColorMode();
% 
% Implements a model for a block design fMRI experiment in which the color
% direction and temporal frequency of stimuli are varied in a block design.
%
% 7/25/16  dhb, bc  Started on this.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods  
    end
    
    properties (Dependent)
    end
    
    % Methods that must be implemented (Abstract in parent class).
    methods (Access=public)             
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
