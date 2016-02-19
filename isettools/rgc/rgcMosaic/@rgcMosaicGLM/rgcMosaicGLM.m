classdef rgcMosaicGLM < rgcMosaic
%% Define an rgcMosaic cell type whose computation is GLM
%
% This class (@rgcMosaicGLM) is a subclass of @rgcMosaic. It is called when
% creating a new rgcMosaic from an inner retina object.  Typically we get
% here from rgcMosaicCreate with the call:
% 
%      mosaicGLM = rgcMosaicGLM(ir, mosaicType);
%  
% Inputs: 
%       os: an isetbio outer segment structure
%       mosaicType: 'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget', 'Small Bistratified' 
% 
% Outputs: the rgcMosaicGLM object; rgcMosaicCreate attaches the
%       rgcMosaicGLM object to an innerRetina object.
% 
% The coupled GLM model implemented here is found in Pillow, Shlens,
% Paninski, Sher, Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% This model incorporates code by Pillow available at 
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
% 
% Example: from rgcGLM.m initiailize:
%        
%       os  = osCreate('identity');
%       innerRetina = irCreate(os,'GLM','name','myRGC'); 
%       innerRetina.mosaicCreate('model','glm','mosaicType','on midget');
% 
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
%         numberTrials;
        
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
