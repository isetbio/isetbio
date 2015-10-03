classdef rgcMosaic < handle
% @rgcMosaic: a class with a parent @rgc object, each @rgcMosaic stores
% properties for one type of RGC cell mosaic, e.g., ON midget or small
% bistratified. The parent @rgc object stores each @rgcMosaic as a cell in
% rgcObj.mosaic such that the properties of an @rgcMosaic object can be
% viewed with 'rgcObj.mosaic{1}'. 
% 
% cellTypeInd = 1: ON parasol
% cellTypeInd = 2: OFF parasol
% cellTypeInd = 3: ON midget
% cellTypeInd = 4: OFF midget
% cellTypeInd = 5: small bistratified
    
% 9/2015 JRG

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = protected, GetAccess = public)
        parent;
        input;
        nameCellType;
        receptiveFieldDiameter1STD;
        spatialRFArray;
        spatialRFcenter;
        spatialRFsurround;
        spatialRFonedim;
        spatialRFcontours;
        spatialRFFill;
        cellCenterLocations;
        temporalImpulseResponseCenterRGB;
        temporalImpulseResponseSurroundRGB;
        linearResponse;

    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcMosaic(rgc, sensor, outersegment, varargin)
            % Initialize the parent class
            %

            % Initialize ourselves
            obj.initialize(rgc, sensor, outersegment, varargin{:});
            
            % % parse the varargin
            % for k = 1:2:numel(varargin)
            %     obj.(varargin{k}) = varargin{k+1};
            % end
        end
        
        % set function, see for details
        function obj = set(obj, param, val, varargin)
            mosaicSet(obj, param, val, varargin);
        end
        
        % get function, see for details
        function val = get(obj, param, varargin)
           val = mosaicGet(obj, param, varargin);
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
%         function obj = compute(obj, sensor, outersegment, varargin)
%             % see for details
%             % obj = mosaicCompute(obj, sensor, outersegment, varargin); 
%         end
%         function plot(obj, sensor)
%             % see for details
%             % mosaicPlot(obj, sensor);
%         end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj, sensor, outersegment, varargin);
    end
    
end
