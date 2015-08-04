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
        
        % set function, see osLinearSet for details
        function obj = set(obj, param, val, varargin)
            osLinearSet(obj, param, val, varargin);
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, param, varargin)
           val = osLinearGet(obj, param, varargin);
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sensor, param, varargin)
            osCompute(obj, sensor, varargin); 
        end
        plotResults(obj, sensor);
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        filterKernel(obj, varargin);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
