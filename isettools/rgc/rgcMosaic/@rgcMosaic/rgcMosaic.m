classdef rgcMosaic < handle
%% rgcMosaic stores properties for one type of RGC cell mosaic
% e.g., ON midget or small bistratified. 
% 
% The RGC class stores multiple mosaics as a cell array,
% as well as a few other basic properties about the inner retina
% position and the stimulus size.
% 
% cellTypeInd = 1: ON parasol
% cellTypeInd = 2: OFF parasol
% cellTypeInd = 3: ON midget
% cellTypeInd = 4: OFF midget
% cellTypeInd = 5: small bistratified
% 
% Properties:
%         cellType: one of the above five strings
%         rfDiameter: the 1 stdev RF diameter of RGCs in micrometers
%         rfDiaMagnitude: the magnitude of the RF at 1 stdev distance
%         cellLocation: the spatial location of the center of the RF
%         sRFcenter: the spatial RF center, where (0,0) is at cellLocation
%         sRFsurround: the spatial RF surround
%         tCenter: the temporal impulse response of the center (in RGB)
%         tSurround: the temporal impulse response of the surround (in RGB)
%         linearResponse: the result of a linear convolution of the sRF and
%           the tRF with the RGB stimulus.
% 
% 9/2015 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
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
    linearResponse;     % Store the linear response after convolution
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
