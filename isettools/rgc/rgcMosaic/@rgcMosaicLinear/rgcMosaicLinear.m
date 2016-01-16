classdef rgcMosaicLinear < rgcMosaic
% @rgcMosaicLinear: a subclass of @rgcMosaic. This function is only called by
% rgcLinear to initiailize a mosaic of the rgc object.
% 
%        rgc.mosaic{ind} = rgcMosaicLinear(cellTypeInd, rgc, scene, sensor, outersegment, varargin{:});
% 
% @rgcLinear: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. The
%l inear model follows the details outlined in
% Chichilnisky & Kalmar (2002), and incorporates other anatomical and
% physiological data from several other sources for parameters like
% receptive field spacing, spatial/temporal linear filters and
% nonlinearities. See comments below for details and references.
%
% Inputs: 
%       scene: an isetbio scene structure
%       sensor: an isetbio sensor structure
%       os: an isetbio outer segment structure
%    Optional but recommended:
%       eyeSide: 'left' or 'right', which eye the retinal patch is from
%       patchRadius: radius of retinal patch in microns
%       patchAngle: polar angle of retinal patch
%     [These inputs determine the size of spatial receptive fields, and are
%       necessary to accurately model physiological responses.]
% 
% Outputs: the rgc object.
% 
% Models found in Chichilnisky & Kalmar, J. Neurosci (2002).
% 
% Example: from rgcLinear.m initiailize:
%        obj.mosaic{cellTypeInd} = rgcMosaicLinear(cellTypeInd, obj, scene, sensor, outersegment, varargin{:});
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
        function obj = rgcMosaicLinear(rgc, cellTypeInd, outersegment, sensor, scene, varargin)
            % Initialize the parent class
            obj = obj@rgcMosaic(rgc, cellTypeInd, outersegment, sensor, scene, varargin{:});
            
        end
        
        % set function, see for details
        function obj = set(obj, varargin)
            mosaicSet(obj, varargin{:});
        end
        
        % get function, see for details
        function val = get(obj, varargin)
           val = mosaicGet(obj, varargin{:});
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
