classdef osIdentity < outerSegment 
% @osIdentity: a subclass of @outerSegment object
% 
% This subclass bypass the temporal filtering of the other outer segment
% subclasses and passes the cone isomerizations without modification. 
% It is intended to be used for stimulus-referred retinal
% ganglion cell models initially.
% 
% identityOS = osIdentity(); 
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
        function obj = compute(obj, sensor)
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
