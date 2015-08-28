classdef osLinear < outerSegment 
% @osLinear: Subclass of @outerSegment object
% 
% This subclass implements the cone linear filters determined by the
% experiments found in Angueyra and Rieke (2013) to convert cone 
% isomerizations (R*) to current (pA).
% 
% The osLinear object calculates the outer segment 
% current by convolving linear filters for the L, M and S cones with the 
% isomerization signal.
% 
% linearOS = osLinear(); 
%
% 7/2015 JRG

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % These are the linear filters generated below via filterKernel.
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
            % see osLinearCompute for details
            obj = osLinearCompute(obj, sensor, varargin); 
        end
        function plot(obj, sensor)
            % see osLinearPlot for details
            osLinearPlot(obj, sensor);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        % Generates the linear filters for L, M and S cones using data from
        % physiology by Angueyra and Rieke (2013)
        filterKernel(obj, varargin);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
