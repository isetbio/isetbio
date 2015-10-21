classdef outerSegment < handle
% The @outerSegment parent class for modeling the responses 
% of the cone outer segment. The subclasses @osLinear and @osBioPhys
% calculate outer segment current responses using either a linear temporal
% filter model or a biophysical difference equations model (both based on
% physiology from Fred Rieke's lab, and both models determined by Fred
% Rieke). The class also allows the addition of noise based on
% physiological measurements of cone responses to a neutral gray background
% (physiology and model also by Fred Rieke). 
%
% See subclasses osLinear and osBioPhys for examples.
% 
% JRG, NC, DHB, 8/2015
% 
    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    % The noiseFlag determines if the noisy version of the current is
    % calculated.
    % The current signal is stored in ConeCurrentSignal in an array the
    % same size as the cone mosaic, and if the noiseFlag is set to 1, the
    % noisy version is stored in ConeCurrentSignalPlusNoise.
    properties (SetAccess = protected)

        noiseFlag;
        ConeCurrentSignal;
        ConeCurrentSignalPlusNoise;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    % Public methods
    methods
        
        function obj = outerSegment(varargin)
            % obj.initialize(varargin{:});
            
            obj.noiseFlag = 0;            
            obj.ConeCurrentSignal = [];
            obj.ConeCurrentSignalPlusNoise = [];
        end
        
        % see osSet in @osLinear and @osBioPhys for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % see osGet in @osLinear and @osBioPhys for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end
        
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % see osLinearCompute, osBioPhysCompute
        compute(obj, sensor, param, varargin);
        % see osLinearPlot, osBioPhysPlot
        plot(obj, sensor);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
