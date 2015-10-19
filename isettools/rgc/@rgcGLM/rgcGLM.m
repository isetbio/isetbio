classdef rgcGLM < rgc 
% @rgcGLM: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. The
% GLM (generalized linear model) follows the details outlined in
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, 
% Nature (2008), and incorporates other anatomical and
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
        function obj = rgcGLM(scene, sensor, outersegment, varargin)
            % Initialize the parent class
            obj = obj@rgc(scene, sensor, outersegment, varargin{:});
            
            % Initialize ourselves by building GLM mosaic objects
            for cellTypeInd = 1:length(obj.mosaic)
                obj.mosaic{cellTypeInd} = rgcMosaicGLM(cellTypeInd, obj, scene, sensor, outersegment, varargin{:});
            end
            
        end
        
        % set function, see superclass method in @rgc for details
        function obj = rgcSet(obj, varargin)
            rgcSet@rgc(obj,varargin{:});
        end
        
        % get function, see superclass method in @rgc for details
        function val = rgcGet(obj, varargin)
           % val = rgcGet(obj, varargin{:});
           val = rgcGet@rgc(obj,varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = rgcCompute(obj, outersegment, varargin)
            obj = rgcCompute@rgc(obj, outersegment, varargin{:});
        end
        function rgcPlot(obj, varargin)
            rgcPlot@rgc(obj, varargin{:});
        end
        function rgcMovie(obj, outersegment, varargin)
            rgcMovie@rgc(obj, outersegment, varargin{:});
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
