classdef osIdentity < outerSegment 
% @osIdentity: a subclass of @outerSegment object
% 
% This subclass implements the cone Identity filters determined by the
% experiments found in Angueyra and Rieke (2013) to convert cone 
% isomerizations (R*) to current (pA).
% 
% The osIdentity object calculates the outer segment 
% current by convolving Identity filters for the L, M and S cones with the 
% isomerization signal.
% 
% IdentityOS = osIdentity(); 
%
% 7/2015 JRG

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % These are the Identity filters generated below via filterKernel.
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = osIdentity(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
            
            % Initialize ourselves
            obj.initialize();
            
            % parse the varargin
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
        end
        
        % set function, see osIdentitySet for details
        function obj = set(obj, param, val, varargin)
            osIdentitySet(obj, param, val, varargin);
        end
        
        % get function, see osIdentityGet for details
        function val = get(obj, param, varargin)
           val = osIdentityGet(obj, param, varargin);
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sensor, param, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin); 
        end
        function plot(obj, sensor)
            % see osPlot for details
            osPlot(obj, sensor);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
