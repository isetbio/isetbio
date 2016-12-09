classdef bipolar < handle
% Define the bipolar cell class.
% 
% The bipolar class allows the simulation of retinal processing from the
% cone outer segment current to the retinal ganglion cell spike response.
% Although we do not yet have a fully validated model of this architecture,
% this represents a first attempt at simulating RGC responses from the cone
% photocurrent. 
% 
% In order to achieve this, the bipolar temporal impulse response (IR)
% is the result of deconvolution of the cone IR from the RGC IR.
%
% The bipolar object also allows for the simulation of nonlinear subunits
% within retinal ganglion cell spatial receptive fields.
% 
% Input: the cone response over time from an outer segment object.
% 
% Output: the bipolar response over time, which can be fed into an inner
% retina object.
% 
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    
    cellLocation;                    % location of bipolar RF center
    cellType;                        % diffuse on or off
    patchSize;                       % size of retinal patch from sensor
    timeStep;                        % time step of simulation from sensor
    filterType;                      % bipolar temporal filter type
    sRFcenter;                       % spatial RF of the center on the receptor grid
    sRFsurround;                     % spatial RF of the surround on the receptor grid
    rectificationCenter              % nonlinear function for center
    rectificationSurround            % nonlinear function for surround
    responseCenter;                  % Store the linear response of the center after convolution
    responseSurround;                % Store the linear response of the surround after convolution

end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
    coneType;
end

% Public methods
methods
    
    % Constructor
    function obj = bipolar(cmosaic, varargin)     
        
        p = inputParser;
        addRequired(p,  'cmosaic');
        addParameter(p, 'cellType', 'offDiffuse', @ischar);
        addParameter(p, 'rectifyType', 1, @isnumeric);
        addParameter(p, 'filterType',  1, @isnumeric);
        addParameter(p, 'cellLocation',  [], @isnumeric);
        addParameter(p, 'ecc',  1, @isnumeric);
        addParameter(p, 'coneType',  -1, @isnumeric);
        
        p.parse(cmosaic, varargin{:});  
        
        % Store the spatial pattern of input cones
        obj.coneType  = cmosaic.pattern;
        
        % Match the time step of the cone mosaic
        os = cmosaic.os;
        obj.patchSize = osGet(os,'patchSize');
        obj.timeStep  = cmosaic.integrationTime;
        
        obj.cellType = p.Results.cellType;
        
        switch p.Results.rectifyType
            case 1
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) zeros(size(x));
            case 2
                obj.rectificationCenter = @(x) x.*(x>0);
                obj.rectificationSurround = @(x) zeros(size(x));                
            case 3
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) x;
            case 4
                obj.rectificationCenter = @(x) x.*(x>0);
                obj.rectificationSurround = @(x) x.*(x<0);
            otherwise
                
                obj.rectificationCenter = @(x) x;
                obj.rectificationSurround = @(x) zeros(size(x));
        end
        
        obj.filterType = p.Results.filterType;
        
        % Build spatial receptive field
        
        % The data for the size of the support is based off of this passage
        % from Dacey, Brainard, Lee, et al., Vision Research, 2000, page
        % 1808 bottom right.
        % (http://www.cns.nyu.edu/~tony/vns/readings/dacey-etal-2000.pdf)
        %
        % Also see Boycott & Wassle, 1991,  (European Journal of
        % Neuroscience), Table 1.
        %
        % Likely the larger RF sizes measured physiological (Dacey et al.)
        % vs anatomically (B&W) reflect spread of signals among cones (via
        % direct gap junctions) and probably more important among cone
        % bipolars (via gap junctions with AII amacrine cells). - Fred
        
        % @JRG included  preferential cone selection to bipolar rfs. These
        % set the basic parameters of the spatial receptive fields. There
        % are some detailed modifications about the type of cone inputs (no
        % S-cones to the on midget or the diffuse systems).
        % 
        % We will also incorporate a function that changes the size of the
        % spread and support as a function of eccentricity.  For now we
        % just put in some placeholder numbers.
        
        % These numbers don't make sense to BW at this time.  We
        % need to write a script showing how big they are with
        % respect to the cone mosaic, and we need to check how they
        % vary with eccentricity.  Comparing with the curves in the
        % cited data would be best.
        
        switch obj.cellType
            case{'onDiffuse','offDiffuse'}
                % Diffuse bipolars that carry parasol signals
                % ecc = 0 mm yields 2x2 cone input to bp
                % ecc = 30 mm yields 5x5 cone input to bp
                
                % BW, screwing around.  Just made spatial spread up here.
                % Support formula extrapolated from data in Dacey ... Lee, 1999 @JRG to insert
                support = max(7,floor(2 + 3/10*(p.Results.ecc))); 
                spread = 1;  % Standard deviation of the Gaussian - will be a function
                rfCenterBig   = fspecial('gaussian',[support, support],spread);     % convolutional for now
                rfSurroundBig = fspecial('gaussian',[support,support], 1.3*spread); % convolutional for now
                
                obj.sRFcenter   = rfCenterBig(:,:);
                obj.sRFsurround = 0.7*rfSurroundBig(:,:);
                
            case {'onSBC'}
                % Small bistratified cells - handle S-cone signals
                
                % Needs to be checked and thought through some more @JRG
                % for this particular cell type.
                support = floor(2 + 3/10*(p.Results.ecc));
                
                spread = 3;  % Standard deviation of the Gaussian - will be a function
                rfCenterBig   = fspecial('gaussian',[support,support],spread); % convolutional for now
                rfSurroundBig = fspecial('gaussian',[support,support],1.5*spread); % convolutional for now
                
                obj.sRFcenter = rfCenterBig(:,:);
                obj.sRFsurround = rfSurroundBig(:,:);
                
            case{'onMidget','offMidget'}
                % Midget bipolars to midget RGCs
                
                % ecc = 0 mm yields 1x1 cone input to bp
                % ecc = 30 mm yields 3x3 cone input to bp
                % Support formula extrapolated from data in Dacey ... Lee, 1999 @JRG to insert

                support = floor(1 + (2/10)*(p.Results.ecc)); 
                spread = 1;
                obj.sRFcenter   = fspecial('gaussian',[support,support], spread); % convolutional for now
                obj.sRFsurround = fspecial('gaussian',[support,support], 1.5*spread); % convolutional for now
             
        end
        if isfield(p.Results,'cellLocation')
            obj.cellLocation = p.Results.cellLocation;
        end
    end
    
    % see bipolarSet for details
    function obj = set(obj, varargin)
        bipolarSet(obj, varargin{:});
    end
    
    % see bipolarGet for details
    function val = get(obj, varargin)
        val = bipolarGet(obj, varargin{:});
    end
    
    % see bipolarCompute for details
    function val = compute(obj, varargin)
        val = bipolarCompute(obj, varargin{:});
    end
    
    function plot(obj, varargin)
        bipolarPlot(obj, varargin{:});
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
end

end