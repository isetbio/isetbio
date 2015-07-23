classdef osLinear < outerSegment 
% @osLinear: Subclass of @outerSegment that implements a linear
% version of the temporalFilter() method
%

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        sConeFilter;
        mConeFilter;
        lConeFilter;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = osLinear(varargin)
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
        function obj = compute(obj, sensor, param, varargin)
            osCompute(obj, sensor, param, varargin); 
        end
        plotResults(obj, sensor);
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        filterKernel(obj);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
