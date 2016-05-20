classdef osLinear < outerSegment 
% Linear subclass of the outersegment object
%
%   os = osLinear();
% 
% Implements isomerizations (R*) to photocurrent (pA) using only cone
% linear temporal filters.  The default values are those determined by the
% Angueyra and Rieke (2013, Nature Neuroscience).
% 
% The current is calculated by convolving separate temporal filters
% for the L, M and S cones with the isomerization time course.
%
% The time base of the temporal filters is established using parameters
% from sensorGet(sensor, 'time interval').
% 
% See also:
%
% JRG Copyright ISETBIO Team, 2015

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
            obj.matchSensor(varargin{:});
            
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
                
        matchSensor(obj, varargin);

        function obj = compute(obj, sensor, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin{:}); 
        end
        
        function plot(obj, sensor, varargin)
            % see osPlot for details
            osPlot(obj, sensor, varargin{:});
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
