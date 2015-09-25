classdef rgcLinear < rgc 
% @rgcLinear: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. 
% This model follows the details of the linear model outlined in
% Chichilnisky & Kalmar (2002), and incorporates other anatomical and
% physiological data from several other sources for parameters like
% receptive field spacing, spatial/temporal linear filters and
% nonlinearities. See comments below for details and references.
%
% 9/2015 JRG

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)
       

    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcLinear(sensor, outersegment, varargin)
            % Initialize the parent class
            obj = obj@rgc(sensor, outersegment, varargin{:});
            
            % Initialize ourselves
            obj.initialize(sensor, outersegment, varargin{:});
            
            % % parse the varargin
            % for k = 1:2:numel(varargin)
            %     obj.(varargin{k}) = varargin{k+1};
            % end
        end
        
        % set function, see for details
        function obj = set(obj, param, val, varargin)
            rgcSet(obj, param, val, varargin);
        end
        
        % get function, see for details
        function val = get(obj, param, varargin)
           val = rgcGet(obj, param, varargin);
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sensor, outersegment, varargin)
            % see for details
            obj = rgcCompute(obj, sensor, outersegment, varargin); 
        end
        function plot(obj, varargin)
            % see for details
            rgcPlot(obj, varargin);
        end
        function movie(obj, outersegment, varargin)
            rgcMovie(obj, outersegment, varargin)
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj, sensor, outersegment, varargin);
    end
    
end
