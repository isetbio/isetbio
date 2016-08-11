classdef osBioPhys < outerSegment 
% BioPhys subclass of the outersegment object
% 
%       os = osBioPhys();
% 
% Converts isomerizations (R*) to outer segment current (pA). The
% difference equation model by Rieke implements a biophysical
% simulation of the phototransduction cascade. If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
%%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
% 
% JRG/HJ/BW, ISETBIO Team, 2016

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
