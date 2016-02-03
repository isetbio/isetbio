classdef rgcGLM < rgc 
% @rgcGLM: a subclass of @rgc. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. 
% This function is typically called by rgcCreate, but may also be called
% as an alternative to that.
% 
%       rgc = rgcGLM(scene, sensor, outersegment, varargin)
% 
% The GLM (generalized linear model) follows the details outlined in
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, 
% Nature (2008), and incorporates other anatomical and
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
% Models found in Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, 
%       Nature (2008).
% 
% This model incorporates code by Pillow available at
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
% 
% Example:
%       rgc1 = rgcGLM(scene, sensor, outersegment);
% 
%       eyeAngle = 180; % degrees
%       eyeRadius = 3; % mm
%       eyeSide = 'right';
%       rgc2 = rgcGLM(scene, absorptions, os, eyeSide, eyeRadius, eyeAngle);
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
        function obj = rgcGLM(params)
            % Initialize the parent class
            obj = obj@rgc(params);
            
            % Initialize ourselves by building GLM mosaic objects
%             for cellTypeInd = 1:5%length(obj.mosaic)
%                 params.cellTypeInd = cellTypeInd;
%                 obj.mosaic{cellTypeInd,1} = rgcMosaicGLM(obj);
%             end
            
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
        function obj = compute(obj, outersegment, varargin)
            obj = rgcCompute(obj, outersegment, varargin{:});
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
