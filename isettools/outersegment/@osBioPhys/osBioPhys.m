classdef osBioPhys < outerSegment 
% @osBioPhys: a subclass of @outerSegment that implements a biophysical
% model of the phototransduction cascade to convert cone isomerizations
% (R*) to current (pA). The model was built by Fred Rieke, and details can
% be found at:
%
% http://isetbio.github.io/isetbio/cones/adaptation%20model%20-%20rieke.pdf
% and 
% https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% adaptedOS = osBioPhys();
% 
% 7/2015 JRG
    properties(Access = private)
        state;   % biophysics parameter state
    end
    
    methods
        function obj = set(obj, varargin)
            % set function, see osBioPhysSet for details
            osSet(obj, varargin{:});
        end
        
        
        function val = get(obj, varargin)
            % get function, see osBioPhysGet for details
            val = osGet(obj, varargin{:});
        end
    end
    
    methods (Access=public)        
        function obj = compute(obj, sensor, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin{:});
        end
        
        function plot(obj, sensor, varargin)
            % see osPlot for details
            osPlot(obj, sensor, varargin{:});
        end
    end
end
