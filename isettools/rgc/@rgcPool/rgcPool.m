classdef rgcPool < rgc 
% @rgcPool: a subclass of @rgc. 
%
% This subclass implements retinal ganglion cell computations
% with the @outerSegment object as input. This model follows the
% details of the linear model outlined in Chichilnisky & Kalmar
% (2002), and incorporates other anatomical and physiological
% data from several other sources for parameters like receptive field
% spacing for spatial pooling. rgcPool does not perform any temporal
% filtering.
%
% This class has no differences from the superclass of rgc.  So,
% really it is a placeholder for some day when we might decide it
% needs its own properties.
%
% See also:  *rgcCreate* for how to initialize the parameters
%
% 9/2015 JRG Copyright ISETBIO Team


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
        function obj = rgcPool(os, params)
            % Initialize the parameters of the rgc parent class
            % The params are checked before we get here in the
            % rgcCreate() function.
            obj = obj@rgc(os, params);

            % We think people should be forced to decide and
            % choose their mosaics rather than it happening like
            % this, behind the scenes.  To discuss. Initialize
            % the specific linear mosaic properties
%             for cellTypeInd = 1:5%length(obj.mosaic)
%                 params.cellTypeInd = cellTypeInd;
%                 obj.mosaic{cellTypeInd,1} = rgcMosaicPool(obj, params);
%             end
        end
        
        % set function, see superclass method in @rgc for details
        function obj = rgcSet(obj, varargin)
            rgcSet@rgc(obj, varargin{:});
        end
        
        % get function, see superclass method in @rgc for details
        function val = rgcGet(obj, varargin)
           val = rgcGet@rgc(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, outersegment, varargin)
            % see superclass method in @rgc for details
            % obj = rgcCompute@rgc(obj, outersegment, varargin{:}); 
            obj = rgcCompute(obj, outersegment, varargin{:}); 
        end
        function rgcPlot(obj, varargin)
            % see superclass method in @rgc for details
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
