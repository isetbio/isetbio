classdef OuterSegment < handle
% OuterSegment Parent class
%
    % Public read/write properties
    properties
%        inputSignal; 
%        time;
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these
    properties (SetAccess = protected)
%         timeConstant;
%         filterKernel;
%         outputSignal;
        noiseflag
        ConeCurrentSignal
        ConeCurrentSignalPlusNoise
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    % Public methods
    methods
        
        function obj = OuterSegment()
            obj.initialize();
        end
        
        % set function, see outersegmentSet for details
        function obj = set(obj, param, val, varargin)
            outersegmentSet(obj, param, val, varargin);
        end
        
        % method declaration 
        % get function, see outersegmentGet for details
        function val = get(obj, param, varargin)
           val = outersegmentGet(obj, param, varargin);
        end
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        temporalFilter(obj, sensor, param, varargin);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        computeFilter(obj);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
