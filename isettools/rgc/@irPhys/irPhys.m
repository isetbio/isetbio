classdef irPhys < ir 
% @irPhys: a subclass of @ir. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. 
% This function is typically called by irCreate, but may also be called
% as an alternative to that.
% 
%       ir = irPhys(scene, sensor, outersegment, varargin)
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
% Outputs: the ir object.
% 
% Models found in Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, 
%       Nature (2008).
% 
% This model incorporates code by Pillow available at
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
% 
% Example:
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
        function obj = irPhys(outersegment, sensor, varargin)
            % Initialize the parent class
             obj = obj@ir(outersegment, varargin{:});
           % obj = [];
            % Initialize ourselves by building GLM mosaic objects
            for cellTypeInd = 1%:length(obj.mosaic)
                obj.mosaic{cellTypeInd} = rgcMosaicPhys(obj, cellTypeInd, outersegment, sensor, varargin{:});
            end
            
        end
        
        % set function, see superclass method in @ir for details
        function obj = irSet(obj, varargin)
            irSet@ir(obj,varargin{:});
        end
        
        % get function, see superclass method in @ir for details
        function val = irGet(obj, varargin)
           % val = irGet(obj, varargin{:});
           val = irGet@ir(obj,varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, outersegment, varargin)
            obj = irCompute(obj, outersegment, varargin{:});
        end
        function irPlot(obj, varargin)
            irPlot@ir(obj, varargin{:});
        end
        function irMovie(obj, outersegment, varargin)
            irMovie@ir(obj, outersegment, varargin{:});
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        % initialize(obj, sensor, outersegment, varargin);
    end
    
end
