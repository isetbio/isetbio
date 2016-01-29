classdef rgc < handle
% The @rgc parent class is a component of isetbio that is used to compute
% the responses of retinal ganglion cells from the isetbio @scene or @outerSegment 
% objects. Subclasses of @rgc include @rgcLinear, @rgcLNP and @rgcGLM and 
% are used to instantiate specific models of RGC processing.
% 
% See isetbio/scripts/xNeedsChecking/rgc/s_rgc.m for a tutorial.
% 
% See Pillow, Jonathan W., et al. "Spatio-temporal correlations and visual
% signalling in a complete neuronal population." Nature 454.7207 (2008)
% and Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries 
% in ON and OFF ganglion cells of primate retina." The Journal of 
% Neuroscience 22.7 (2002).
% 
% Properties:
%         name:  animal, RGC; example: 'macaque RGC'
%         input: RGB stim or outersegment
%         temporalEquivEcc: calculated from retinal position, see retinalLocationToTEE        
%         mosaic: a cell array, where each cell is an rgcMosaic object, 
%               which is a subclass of the rgc object.
% 
% Methods: intialize, set, get, compute, plot, movie
%       see individual .m files for details.
%        
% 
% (c) isetbio
% 
% 9/2015 JRG
% 
% 

    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        name; 
        input;
        eyeSide;
        eyeRadius;
        eyeAngle;
        temporalEquivEcc;        
        mosaic;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    % Public methods
    methods
        
        function obj = rgc(outersegment, varargin)
            obj.initialize(outersegment, varargin{:});
        end
        
        % set function, see rgcSet
        function obj = set(obj, varargin)
            obj = rgcSet(obj, varargin{:});
        end
        
        % get function, see rgcGet
        function val = get(obj, varargin)
           val = rgcGet(obj, varargin{:});
        end
        
        % compute function, see rgcCompute
        function obj = compute(obj, outerSegment, varargin)
            obj = rgcCompute(obj,  outerSegment, varargin{:});
        end
        
        % plot function, see rgcPlot
        function plot(obj, varargin)
            rgcPlot(obj, varargin{:});
        end
        
        % movie function, see rgcMovie
        function movie(obj, outersegment, varargin)
            rgcMovie(obj, outersegment, varargin{:});
        end
        
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    
    methods (Abstract, Access=public)
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

        spConvolve
        fullConvolve
       
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj, sensor, outersegment, varargin);
    end
    
end

    
    