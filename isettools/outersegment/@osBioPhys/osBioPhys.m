classdef osBioPhys < outerSegment 
% @NonLinearOuterSegment: Subclass of @OuterSegment that implements a nonlinear
% version of the temporalFilter() method
%

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = osBioPhys(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
            
            % Initialize ourselves
            obj.initialize();
            
            % parse the varargin
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)        
        function obj = compute(obj, sensor, varargin)
            osCompute(obj, sensor, varargin);
        end
        plotResults(obj, sensor);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
