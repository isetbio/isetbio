classdef bipolarMosaic < cellMosaic
%BIPOLARMOSAIC - Create a bipolar mosaic object
%
% The bipolar mosaic class is a subclass of cellMosaic. It is used to
% simulate the processing from the cone outer segment current to the
% bipolar cell current output. 
%
% To create the bipolar mosaic, use
% 
%   bp = bipolarMosaic(coneMosaic, 'PARAM1', val1, 'PARAM2', val2,...) 
%
% Optional parameter name/value pairs are listed below.
% 
% The bipolar object also allows for the simulation of nonlinear subunits
% within retinal ganglion cell spatial receptive fields.
% 
% Input (Required):
%  coneMosaic: a cone mosaic object that includes the photocurrent response
% 
% Input (Parameter-Value pairs):
%     cellLocation;                    % location of bipolar RF center
%     cellType;                        % there are five different cell types
%     patchSize;                       % size of retinal patch from sensor
%     timeStep;                        % time step of simulation from sensor
%     filterType;                      % bipolar temporal filter type
%     sRFcenter;                       % spatial RF of the center on the receptor grid
%     sRFsurround;                     % spatial RF of the surround on the receptor grid
%     rectificationCenter              % nonlinear function for center
%     rectificationSurround            % nonlinear function for surround
%     responseCenter;                  % Store the linear response of the center after convolution
%     responseSurround;                % Store the linear response of the surround after convolution  
% 
%  ISETBIO wiki: <a href="matlab:
%  web('https://github.com/isetbio/isetbio/wiki/bipolar','-browser')">bipolar</a>.
%  
% Scientific notes and references
%  We do not have a validated model of the bipolar temporal impulse
%  response.  To simulate it, we assume that the RGC responses from the
%  cone photocurrent as measured by EJ is correct.  The bipolar temporal
%  impulse response (tIR) is the response necessary to combine with the
%  cone temporal tIR to equal the RGC tIR.
%

% JRG/BW (c) isetbio team, 2016

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    % TODO: We should specify an amplitude for the center and surround.
    % Perhaps we should specify parameters of the receptive fields beyond
    % what is in cellMosaic.
    
    % RECTIFICATIONCENTER nonlinear function for center
    rectificationCenter              
    
    % RECTIFICATIONSURROUND nonlinear function for surround
    rectificationSurround            
    
    % RESPONSECENTER Store the linear response of the center after convolution
    responseCenter;                  
    
    % RESPONSESURROUND Store the linear response of the surround after convolution
    responseSurround;

end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
    % CONETYPE  
    % on diffuse, off diffuse and on midget bipolars get no S cone input
    % off midget bipolars get L,M,S cone input to center
    % on sbc bipolars get only S cone center and only L+M surround
    coneType;
end

properties (Access = public)
end

% Public methods
methods
    
    % Constructor
    function obj = bipolarMosaic(cmosaic, varargin)     
        % Initialize the bipolar class
        %
        % Example:
        %   bp = bipolarMosaic(cMosaic,'cellType','ondiffuse');
        %
        
        p = inputParser;
        p.KeepUnmatched;    % Keeps spread and stride 
        p.addRequired('cmosaic',@(x)(isequal(class(x),'coneMosaic')));
        
        p.addParameter('parent',[], @(x)(isequal(class(x),'bipolarLayer')));
        p.addParameter('cellType', 'offdiffuse', @(x)(ismember(ieParamFormat(x),obj.validCellTypes)));
        p.addParameter('rectifyType', 1, @isnumeric);
        p.addParameter('filterType',  1, @isnumeric);
        p.addParameter('ecc',  0, @isnumeric);
        p.addParameter('spread',  1, @isnumeric);
        p.addParameter('stride',  1, @isnumeric);

        p.parse(cmosaic, varargin{:});  
        
        % The layer object that this is part of.
        obj.parent    = p.Results.parent;
        obj.input     = p.Results.cmosaic;
        
        % Store the spatial pattern of input cones
        obj.coneType  = cmosaic.pattern;
        
        % This might be a mistake, but we store the size and time step of
        % the cone mosaic in this mosaic for easy accessibility.  The
        % cMosaic itself is stored in the layer object that contains this
        % mosaic.  Maybe these should just be private variables?
        obj.patchSize = cmosaic.size; 
        obj.timeStep  = cmosaic.integrationTime;
        
        obj.cellType = ieParamFormat(p.Results.cellType);
        
        % Set the rectification operation
        switch p.Results.rectifyType
            case 1
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) zeros(size(x));
            case 2
                obj.rectificationCenter = @(x) x.*(x>0);
                obj.rectificationSurround = @(x) zeros(size(x));                    
            case 3
                obj.rectificationCenter = @(x) x.*(x<0);
                obj.rectificationSurround = @(x) zeros(size(x));  
            case 4
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) x;
            case 5
                obj.rectificationCenter = @(x) x.*(x>0);
                obj.rectificationSurround = @(x) x.*(x<0);
            otherwise
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) zeros(size(x));
        end
        
        obj.filterType = p.Results.filterType;
        
        % Build spatial receptive fields
        obj.spatialRFInit('conemosaic',cmosaic,...
            'ecc',   p.Results.ecc,...
            'stride',p.Results.stride, ...
            'spread',p.Results.spread);
        
    end
    
    function window(obj)
        % Tip: Retrieve guidata using
        %    gui = guidata(obj.figureHandle);
        obj.fig = bipolarWindow(obj);
    end
    
    function bipolarsPerMicron = cellsPerDistance(obj,varargin)
        p = inputParser;
        p.addRequired('obj',@(x)(isequal(class(obj),'bipolarMosaic')));
        p.addParameter('units','m',@ischar);
        p.parse(obj,varargin{:});
        
        patchSizeUM = obj.Parent.size;   % In microns
        
        % The bipolar mosaics at this point are all the same row/col count.
        % But they may not be in the future.  So, what do we do about that?
        bpRowCol = size(obj.Parent.input.mosaic{1}.cellLocation);    
        
        % Converts a distance in microns to a number of bipolars per micron
        bipolarsPerMicron = bpRowCol(1:2) ./ patchSizeUM;   % cells/micron
    end
    
end

properties (Constant)
    % VALIDCELLTYPES 
    %
    % Cell array of strings containing valid values for the cell type.
    % diffuse and parasol are synonyms.  Not sure we should have them both.
    % And, possibly we should have smallbistratified. (BW)
    validCellTypes = {'ondiffuse','offdiffuse','onparasol','offparasol','onmidget','offmidget','onsbc'};
end

% Methods that must only be implemented (Abstract in parent class).

% Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
end

end