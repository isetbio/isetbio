classdef outersegment < handle
    %Class outersegment 
    %   Detailed explanation goes here
    %   Example implementation code here
    %   James Golden, 6/22/15
    
    properties (Access = public)
    end
    
    properties (SetAccess = private, GetAccess = public)
        noiseflag
        ConeCurrentSignal
        ConeCurrentSignalPlusNoise
    end
    
    methods (Access = public)
        
        % constructor
        function obj = outersegment(noiseFlag, varargin)
            % initialize properties
            obj.noiseflag = noiseFlag;
        end
       
        % method declaration
        % set function, see outersegmentSet for details
        function obj = set(obj, param, val, varargin)
            outersegmentSet(obj, param, val, varargin);
        end
        
        % method declaration 
        % get function, see outersegmentGet for details
        function val = get(obj, param, varargin)
           val = outersegmentGet(obj, param, varargin);
        end
        
        % method declaration
        % define as above
        obj = outersegmentCompute(obj, param, sensor, varargin);
    end
     
end

