classdef rgcLNP < rgc 
% @rgcLNP: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. 
% This function is typically called by rgcCreate, but may also be called
% as an alternative to that.
% 
%       rgc = rgcLNP(scene, sensor, outersegment, varargin)
% 
% @rgcLNP: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. The
% LNP (linear-nonlinear-Poisson) model follows the details outlined in
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
% Models found in Chichilnisky & Kalmar, J. Neurosci (2002) and Pillow, Shlens, 
%   Paninski, Sher, Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% This model incorporates code by Pillow available at
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
% 
% Example:
%       rgc1 = rgcLNP(scene, sensor, outersegment);
% 
%       eyeAngle = 180; % degrees
%       eyeRadius = 3; % mm
%       eyeSide = 'right';
%       rgc2 = rgcLNP(scene, absorptions, os, eyeSide, eyeRadius, eyeAngle);
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
        function obj = compute(obj, outersegment, varargin)
            % see for details 
            obj = rgcCompute(obj,  outersegment, varargin{:}); 
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
