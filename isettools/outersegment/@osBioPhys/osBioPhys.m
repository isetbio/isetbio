classdef osBioPhys < outerSegment 
% @osBioPhys: a subclass of @outerSegment that implements a biophysical model
% of the phototransduction cascade to convert cone isomerizations (R*) to 
% current (pA). The model was built by Fred Rieke, and details can be found
% at:
%
% http://isetbio.github.io/isetbio/cones/adaptation%20model%20-%20rieke.pdf
% and 
% https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% adaptedOS = osBioPhys();
% 
% 7/2015 JRG


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
        end
        
        % set function, see osBioPhysSet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osBioPhysGet for details
        function val = get(obj,  varargin)
           val = osGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)        
        function obj = compute(obj, sensor, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin{:});
        end
        function plot(obj, sensor)
            % see osPlot for details
            osPlot(obj, sensor);
        end
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        
    end
    
end
