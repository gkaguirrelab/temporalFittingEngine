classdef tfe < handle
% The @tfe parent class for modeling temporal responses.  This implements
% methods that are model independent and defines the API for subclasses
% that implement particular models.
%
% Methods implemented in this class:
%   fitResponse - Use a model to fit a temporal response
%   fitError - Compute error between model fit and response
%   applyKernal - Apply a convolution kernal to a response
%   concatenatePackets - Take multiple data packets and concatenate into one
%     long packet. 
%   resampleTimebase - Change temporal sampling of a signal
%   plot - Plot response and fit
%
% Methods that must be implemented in subclasses
%   computeResponse - Compute model response
%   isPacket - Verify that packet has right format for model
%   defaultParams - Return a default set of parameters
%   paramsToVec - Convert parameters structure to vector format
%   vecToParams - Conver parameters vector to structure format
%   lockMatrix - Construct parameter locking matrix for the model
%   paramPrint - Print useful things about the parameters
%
% The class constructor can take optional key/value pairs as follows.
%  'verbosity' - string (default 'none').  Verbsoity level for method
%      diagnostic printout
%      'none' - No diagnostic printout
%      'high' - Print everything we can think might be useful.

% 6/26/16  dhb  Started in on this
% 9/21/16  gka  Massive restructuring

    % Public read/write properties
    properties              
    end
    
    % Dependent properties, computed from other parameters
    properties
    end
        
    % Public, read-only properties.  These can be set by methods of the
    % parent class (that is, this class) but not by methods of subclasses.
    properties (SetAccess = private, GetAccess = public)
        % Verbosity level, set at constructor time.
        verbosity
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
    % explicitly decare them.
    %
    % We do stick the constructor method here.
    methods (Access=public)
        % Constructor
        function obj = tfe(varargin)
                     
            %% Parse vargin for options passed here
            p = inputParser;
            p.addParameter('verbosity','none',@ischar);
            p.parse(varargin{:});
            obj.verbosity = p.Results.verbosity;
        end
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % Default parameters and bounds in structure form
        [params,paramsLb,paramsUb] = defaultParams(obj,varargin);
        
        % Print parameters
        paramPrint(obj,params,varagin);
        
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
        
        % Test packets for validity for the model subclass
        packetValidity = isPacket(obj,thePacket)
                
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
