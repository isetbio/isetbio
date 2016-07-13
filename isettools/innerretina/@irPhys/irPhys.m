classdef irPhys < ir 
% @irPhys: a subclass of @ir. This subclass implements retinal
% ganglion cell computations with the @outerSegment object as input. 
% This function is typically called by irCreate, but may also be called
% as an alternative to that.
% 
%       ir = irPhys(outersegment, params)
% 
% The GLM (generalized linear model) follows the details outlined in
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, 
% Nature (2008). irPhys generates a mosaic of RGCs with parameters from GLM
% fits from an experiment in the Chichilnisky lab.
%
% Inputs: 
%       os: an isetbio outer segment structure
%    Optional but recommended:
%       eyeSide: 'left' or 'right', which eye the retinal patch is from
%       patchRadius: radius of retinal patch in microns
%       patchAngle: polar angle of retinal patch
%     [For the irPhys object, these parameters do not have an affect on the
%       properties of the receptive field s, because the GLM parameters
%       are loaded from .mat files. The user should assign the properties
%       according to the values measured for the retinal tissue as recorded
%       in the experiment.]
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
        function obj = irPhys(outersegment, varargin)
            % Initialize the parent class
            obj = obj@ir(outersegment, varargin{:});
            
            % Initialize ourselves by building rgcPhys mosaic objects
            obj.mosaic{1} = rgcPhys(obj, varargin{:});
            
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
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
