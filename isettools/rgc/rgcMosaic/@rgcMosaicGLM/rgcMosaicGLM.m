classdef rgcMosaicGLM < rgcMosaic
% rgcMosaic cell type with a GLM (coupled-nonlinear)computational model
%
% The coupled GLM model is published in Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky & Simoncelli, Nature (2008).% The computational model
% implemented here relies on code by
% <http://pillowlab.princeton.edu/code_GLM.html Pillow>, which is
% distributed under the GNU General Public License.
%
% rgcMosaicGLM is a subclass of rgcMosaic. It is called when creating a new
% rgcMosaic from an inner retina object.  Typically we get here from the
% inner retina object via a call
%
%   ir.mosaicCreate( ... )
% 
% See also:
%
% Example:
%   os = osCreate('identity');        % A pass through from the stimulus
%   ir = irCreate(os,'name','myRGC'); % An inner retina container
%   ir.mosaicCreate('model','glm','mosaicType','on midget'); % This  mosaic
%
% 9/2015 JRG

%% Properties 
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)

        generatorFunction;
        nlResponse;
        numberTrials = 1;
        
        postSpikeFilter;
        couplingFilter;
        couplingMatrix;
        
        spikeResponse;
%         rasterResponse;
%         psthResponse;

    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
  
%% Methods
    % Public methods
    methods
        
        % Constructor
        function obj = rgcMosaicGLM(rgc, mosaicType)
            
            % Initialize for the cell type            
            obj = obj@rgcMosaic(rgc, mosaicType);

            % Initialize for the computationa type, which in this case is
            % always GLM
            obj.initialize(rgc);
%             obj.rgcGLM(rgc);
            
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
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj, varargin);
    end
    
end
