classdef rgcMosaicLNP < rgcMosaic
% @rgcMosaicLNP: a subclass of @rgcMosaic. This function is only called by
% rgcLNP to initiailize a mosaic of the rgc object.
% 
%        rgc.mosaic{ind} = rgcMosaicLNP(cellTypeInd, rgc, scene, sensor, outersegment, varargin{:});
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
% Example: from rgcLNP.m initiailize:
%        obj.mosaic{cellTypeInd} = rgcMosaicLNP(cellTypeInd, obj, scene, sensor, outersegment, varargin{:});
% 
% 9/2015 JRG


    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)

        generatorFunction;
        nlResponse;
        spikeResponse;

        rasterResponse;
        psthResponse;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcMosaicLNP(rgc, sensor, outersegment, varargin)
            % Initialize the parent class            
            obj = obj@rgcMosaic(rgc, sensor, outersegment, varargin{:});

            % Initialize ourselves
            obj.initialize(rgc, sensor, outersegment, varargin{:});
            
            % % parse the varargin
            % for k = 1:2:numel(varargin)
            %     obj.(varargin{k}) = varargin{k+1};
            % end
        end
        
        % set function, see for details
        function obj = set(obj, varargin)
            % obj = set@rgcMosaic(obj, varargin);
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
