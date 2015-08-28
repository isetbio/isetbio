classdef outerSegment < handle
% OuterSegment parent class
%
% Outersegment parent class for modeling the responses 
% of the cone outer segment.
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
        
        function obj = outerSegment()
            obj.initialize();
        end
        
        % set function, see outersegmentSet for details
        % see osLinearSet, osBioPhysSet
        function obj = set(obj, param, val, varargin)
            outersegmentSet(obj, param, val, varargin);
        end
        
        % get function, see outersegmentGet for details
        % see osLinearGet, osBioPhysGet
        function val = get(obj, param, varargin)
           val = outersegmentGet(obj, param, varargin);
        end
        
        function [osCurrent, obj] = adapt(sensor, typeAdapt)
            [osCurrent, obj] = coneAdaptAlt(sensor, typeAdapt);
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
        % called by osLinearCompute
        filterKernel(obj);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
