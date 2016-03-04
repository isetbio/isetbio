classdef rgcMosaic < handle
%% Define the rgcMosaic parent class to store a mosaic of cells with a type and model
%
% This class is called when creating a new rgcMosaic from an inner
% retina object.  Typically we get here from rgcMosaicCreate with the call:
% 
%      mosaicLinear = rgcMosaicLinear(ir, mosaicType);
%      mosaicGLM    = rgcMosaicGLM(ir, mosaicType);
%  
% Inputs: 
%       os: an isetbio outer segment structure
%       mosaicType: 'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget', 'Small Bistratified' 
% 
% Outputs: the rgcMosaic object; rgcMosaicCreate attaches the
%       rgcMosaic object to an innerRetina object.
% 
% The center and surround spatial receptive fields and the temporal impulse
% responses are initialized in @rgcMosaic/initialize. The rgcMosaic
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

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    
    % The type of computational model for the RGC spikes
    % model; 
    
    cellType;           % Possible types are ...
    rfDiameter;         % receptive field center diameter
    % We should estimate the rf center sigma
    % rfDiaMagnitude;
    
    % Cell array cellLocation{i}{j} = [x,y] position
    % (microns)
    cellLocation;
    sRFcenter;          % spatial RF of the center on the receptor grid
    sRFsurround;        % spatial RF of the surround
    tCenter;            % temporal impulse response of the center
    tSurround;          %    and of the surround (1 ms timing by default)
    rfDiaMagnitude;     % for making movies of response
    responseLinear;     % Store the linear response after convolution
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

% Public methods
methods
    
    % Constructor
    function obj = rgcMosaic(rgc, mosaicInd)
        
        % Maybe move all the initialization into here?
        % Initialize ourselves
        obj.initialize(rgc, mosaicInd);
        
    end
    
    % set function, see mosaicSet for details
    function obj = set(obj, param, val, varargin)
        mosaicSet(obj, param, val, varargin{:});
    end
    
    % get function, see mosaicGet for details
    function val = get(obj, param, varargin)
        val = mosaicGet(obj, param, varargin{:});
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
