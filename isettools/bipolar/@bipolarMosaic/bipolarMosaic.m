classdef bipolarMosaic < handle
%BIPOLARMOSAIC - Create a bipolar mosaic object
%
% The bipolar mosaic class allows the simulation of retinal processing from
% the cone outer segment current to the retinal ganglion cell spike
% response.
% 
%    bp = bipolarMosaic(cMosaic, 'PARAM1', val1, 'PARAM2', val2,...) 
%
% creates the bipolar object. Optional parameter name/value pairs are
% listed below.
% 
% We do not yet have a validated model of the bipolar temporal impulse
% response.  To simulate it, we assume that the RGC responses from the cone
% photocurrent as measured by EJ is correct.  The bipolar temporal impulse
% response (tIR) is the response necessary to combine with the cone
% temporal tIR to equal the RGC tIR.
%
% The bipolar object also allows for the simulation of nonlinear subunits
% within retinal ganglion cell spatial receptive fields.
% 
% Input: the cone photocurrent response over time from a cone mosaic.
% 
% Output: the bipolar voltage response over time. This is fed into an inner
% retina object. (See s_vaRGC.m for an example.  More will appear).
% 
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
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    % We need an amplitude for the center and surround
    
    % CELLTYPE diffuse on or off
    cellType;                        
    
    %CELLLOCATION location of bipolar RF centers w.r.t. the input samples
    cellLocation;
    
    % PATCHSIZE size of retinal patch w.r.t. the cone mosaic
    patchSize;                       
    
    % TIMESTEP time step of simulation from original cone mosaic
    timeStep;       
    
    % FILTERTYPE bipolar temporal filter type
    filterType; 
    
    % SRFCENTER spatial RF of the center on the input samples
    % Represented w.r.t the input sampling grid
    sRFcenter = [];                       
    
    % SRFSURROUND spatial RF of the surround on the input samples
    % Represented w.r.t the input sampling grid
    sRFsurround = [];                   
    
    % RECTIFICATIONCENTER nonlinear function for center
    rectificationCenter              
    
    % RECTIFICATIONSURROUND nonlinear function for surround
    rectificationSurround            
    
    % RESPONSECENTER Store the linear response of the center after convolution
    responseCenter;                  
    
    % RESPONSESURROUND Store the linear response of the surround after convolution
    responseSurround;
    
    % cone mosaic input
    input;
    
    % parent - the bipolarLayer containing this mosaic
    parent;

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
    %FIGUREHANDLE When we open the figure for the mosaic, we store the handle here
    figureHandle;
end

% Public methods
methods
    
    % Constructor
    function obj = bipolarMosaic(cmosaic, varargin)     
        % Initialize the bipolar class
        %
        %   bp = bipolar(cMosaic,'cellType','ondiffuse');
        %
        
        p = inputParser;
        addRequired(p,  'cmosaic');
        
        addParameter(p, 'parent',[], @(x)(isequal(class(x),'bipolarLayer')));
        addParameter(p, 'cellType', 'offdiffuse', @(x)(ismember(strrep(lower(x),' ',''),obj.validCellTypes)));
        addParameter(p, 'rectifyType', 1, @isnumeric);
        addParameter(p, 'filterType',  1, @isnumeric);
        addParameter(p, 'cellLocation',  [], @isnumeric);
        addParameter(p, 'ecc',  1, @isnumeric);
        addParameter(p, 'coneType',  -1, @isnumeric);

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
        
        
        obj.cellType = strrep(lower(p.Results.cellType),' ','');
        
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
        
        % Build spatial receptive field
        obj.spatialRFInit('ecc',p.Results.ecc,'conemosaic',cmosaic);
        
    end
    
    function window(obj)
        % Tip: Retrieve guidata using
        %    gui = guidata(obj.figureHandle);
        obj.figureHandle = bipolarWindow(obj);
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
    % VALIDCELLTYPES Cell array of strings containing valid values for the
    % cell type.  diffuse and parasol are synonyms.  Not sure we should
    % have them both.  And, possibly we should have smallbistratified. (BW)
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