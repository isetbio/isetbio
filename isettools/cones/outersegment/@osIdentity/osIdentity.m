classdef osIdentity < outerSegment 
% @osIdentity: A subclass of @outerSegment object
% 
% Syntax:
%	identityOS = osIdentity(); 
%
% Description:
%    This subclass bypass the temporal filtering of the other outer segment
%    subclasses and passes the cone isomerizations without modification. It
%    is intended to be used for stimulus-referred retinal ganglion cell
%    models initially.
%
%    Examples are contained in the code. To access, type 'edit
%    osIdentity.m' into the Command Window.
%
% Inputs:
%    None required.
%
% Outputs:
%    The created outerSegment object
%
% Optional key/value pairs:
%    None.
%

% History:
%    07/xx/15  JRG  Created
%    02/13/18  jnm  Formatting

% Examples:
%{
	identityOS = osIdentity(); 
%}

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        %photonRate - The photon rate
        photonRate
    end
    
    methods    
        % Constructor
        function obj = osIdentity(varargin)
            % Initialize the parent class
            %
            % Syntax:
            %   obj = osIdentity([varargin])
            %
            % Description:
            %    Initialize the parent function.
            %
            % Inputs:
            %    varargin - (Optional) Any additional parameters provided.
            %
            % Outputs:
            %    obj      - The created outersegment identity object.
            %
            % Optional key/value pairs:
            %    None.
            %
            
            % Initialize the parent class
            obj = obj@outerSegment();
            
            % Initialize ourselves
            obj.photonRate = [];
        end
        
        % set function, see osIdentitySet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osIdentityGet for details
        function val = get(obj, varargin)
            val = osGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sceneRGB)
            obj = osCompute(obj, sceneRGB); 
        end
        
        function plot(obj, sensor)
            osPlot(obj, sensor);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        
    end
    
end
