classdef tfeQuadraticColorModel < tfe
% tfeQuadraticColorModel
%
%   tfe = tfeQuadraticColorModel();
% 
% Implements a model that is quadratic in the color contrast of the
% stimulus.
%
% 6/26/16  dhb  Started in on this.

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
