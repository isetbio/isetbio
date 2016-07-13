classdef rgcPhys % < rgcMosaic
% rgcPhys is used for unit testing a GLM computational model for RGC
% response that loads parameters from a physiology experiment in the
% Chichilnisky Lab.
%
% The coupled GLM model is published in Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% The computational model implemented here relies on code by
% <http://pillowlab.princeton.edu/code_GLM.html Pillow>, which is
% distributed under the GNU General Public License.
%
% rgcPhys is not a subclass of rgcMosaic, but is similar to rgcGLM in many
% respects. It is called when creating a new Phys model for an
% inner retina object.  Typically we get here from the inner retina object
% via a call
%
%   ir.mosaicCreate('model','phys','type','your type goes here')
% 
% See also: irPhys, v_rgcExternal, t_rgcNaturalScenesFigure2
%
% Example:
%   os = osCreate('identity');        % A pass through from the stimulus
%   ir = irCreate(os,'name','myRGC'); % An inner retina container
%   ir.mosaicCreate('model','phys'); % This  mosaic
%
% 9/2015 JRG    
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = protected, GetAccess = public)
        experimentID;
        cellType;
        cellID;
        stimulusFit;
        stimulusTest;
        rfDiameter;
        rfDiaMagnitude;
        cellLocation;
        sRFcenter;
        sRFsurround;
        tCenter;
        tSurround;
        responseLinear;
        
        generatorFunction;
        nlResponse;
        numberTrials;
        spikeResponse;
        
        postSpikeFilter;
        couplingFilter;
        couplingMatrix;
        tonicDrive;
        responseRaster;
        responsePsth;
        responseVoltage;
        responseSpikes;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcPhys(rgc, varargin)
            obj = obj.initialize(rgc, varargin{:});            
        end
        
        % set function, see for details
        function obj = set(obj, varargin)
            % obj = set@rgcMosaic(obj, varargin);
            obj = mosaicSet(obj, varargin{:});
        end
        
        % get function, see for details
        function val = get(obj, varargin)
           val = mosaicGet(obj, varargin{:});
        end
        
        [stimX, stimY, offset] = obj.stimPositions(xcell,ycell)
      
        mosaicGLM = obj.glmLoad(cellType);
        
        obj = obj.mosaicLoadExperimental(mosaicGLM, cellType, varargin);
        
        obj = obj.mosaicLoadAverage(mosaicGLM, cellType, varargin);
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        obj = initialize(obj, varargin);
    end
    
end
