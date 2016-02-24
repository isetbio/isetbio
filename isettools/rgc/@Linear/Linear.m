classdef rgcMosaicLinear < rgcMosaic
%% Define an rgcMosaic cell type whose computation is linear
%
% This class (@rgcMosaicLinear) is a subclass of @rgcMosaic. It is called when
% creating a new rgcMosaic from an inner retina object.  Typically we get
% here from rgcMosaicCreate with the call:
% 
%      mosaicLinear = rgcMosaicLinear(ir, mosaicType);
%  
% Inputs: 
%       os: an isetbio outer segment structure
%       mosaicType: 'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget', 'Small Bistratified' 
% 
% Outputs: the rgcMosaicLinear object; rgcMosaicCreate attaches the
%       rgcMosaicLinear object to an innerRetina object.
% 
% The center and surround spatial receptive fields and the temporal impulse
% responses are initialized in @rgcMosaic/initialize. The rgcMosaicLinear
% class is only initialized by the parent class rgcMosaic; it does not have
% its own initialize function. The model implemented here is described in
% Chichilnisky & Kalmar, J. Neurosci (2002).
% 
% Example: from t_rgc.m:
%        
%       os  = osCreate('identity');
%       innerRetina = irCreate(os,'linear','name','myRGC'); 
%       innerRetina.mosaicCreate('model','linear','mosaicType','on midget');
% 
% (c) isetbio team
% 9/2015 JRG
% 
%% Properties
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
        function obj = rgcMosaicLinear(rgc, mosaicInd)
            % Initialize the parent class
            obj = obj@rgcMosaic(rgc, mosaicInd);
            
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
%% Methods
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
