classdef osDisplayRGB < outerSegment 
% displayRGB subclass of the outersegment object
%
%   os = osDisplayRGB();
%
% Bypasses the temporal filtering of the other outer segment subclasses and
% passes the cone isomerizations without modification. It is intended to be
% used for stimulus-referred retinal ganglion cell models initially.
% 
% displayRGBOS = osDisplayRGB(); 
% 
% See also subclasses:
%       osLinear.m, osBioPhys.m
%
% JRG/HJ/BW, ISETBIO Team, 2016

    % Public properties.
    properties (SetAccess = public, GetAccess = public)

    end
    

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        rgbData
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = osDisplayRGB(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
            
            % Initialize ourselves
            obj.rgbData = [];
        end
        
        % set function, see osIdentitySet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osIdentityGet for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sceneRGB)
            % see osCompute for details
            obj = osCompute(obj, sceneRGB); 
        end
        
        function plot(obj, sensor)
            % see osPlot for details
            osPlot(obj, sensor);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        
    end
    
end
