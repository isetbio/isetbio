classdef LinearOuterSegment < OuterSegment 
% @LinearOuterSegment: Subclass of @OuterSegment that implements a linear
% version of the temporalFilter() method
%

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        sConeFilter
        mConeFilter
        lConeFilter
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = LinearOuterSegment(varargin)
            % Initialize the parent class
            obj = obj@OuterSegment();
            
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
        temporalFilter(obj, sensor, param, varargin);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
