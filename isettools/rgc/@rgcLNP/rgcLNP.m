classdef rgcLNP < rgc 
% @rgcLNP: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. The
% LNP (linear-nonlinear-Poisson) model follows the details outlined in
% Chichilnisky & Kalmar (2002), and incorporates other anatomical and
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
        function obj = rgcLNP(scene, sensor, outersegment, varargin)
            % Initialize the parent class
            obj = obj@rgc(scene, sensor, outersegment, varargin{:});
            
            % Initialize ourselves by building LNP mosaic objects
            for cellTypeInd = 1:length(obj.mosaic)
                obj.mosaic{cellTypeInd} = rgcMosaicLNP(cellTypeInd, obj, scene, sensor, outersegment, varargin{:});
            end
        end
        
        % set function, see for details
        function obj = rgcSet(obj, param, val, varargin)
            rgcSet@rgc(obj, param, val, varargin{:});
        end
        
        % get function, see for details
        function val = rgcGet(obj, param, varargin)
           val = rgcGet@rgc(obj, param, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = rgcCompute(obj, outersegment, varargin)
            % see for details
            obj = rgcCompute@rgc(obj,  outersegment, varargin{:}); 
        end
        function rgcPlot(obj, varargin)
            % see for details
            rgcPlot@rgc(obj, varargin{:});
        end
        function rgcMovie(obj, outersegment, varargin)
            rgcMovie@rgc(obj, outersegment, varargin{:})
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
