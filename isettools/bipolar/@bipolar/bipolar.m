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
% ought to be the result of deconvolution of the cone IR from the RGC IR.
% However, we also want to simulate the nonlinear outer segment model, so
% we approximate this bipolar temporal filter with a weighted
% differentiator, such that responseBipolar(t) = w1*responseCone(t) +
% w2*responseCone'(t), where w1 is close to zero and w2 is a large negative
% number.

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
    coneMosaic;
end

% Public methods
methods
    
    % Constructor
    function obj = bipolar(os, varargin)     
        
        p = inputParser;
        addRequired(p, 'os');
        addParameter(p, 'cellType', 'offDiffuse', @ischar);
        addParameter(p, 'rectifyType', 1, @isnumeric);
        addParameter(p, 'filterType',  1, @isnumeric);
        addParameter(p, 'cellLocation',  [], @isnumeric);
        addParameter(p, 'ecc',  1, @isnumeric);
        addParameter(p, 'coneMosaic',  ones(size(os.coneCurrentSignal)), @isnumeric);
        
        p.parse(os, varargin{:});  
        
        obj.patchSize = osGet(os,'patchSize');
        obj.timeStep = osGet(os,'timeStep');
        
        obj.cellType = p.Results.cellType;
        obj.coneMosaic = p.Results.coneMosaic;
        
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
        
        % This plot is based off of this passage from Dacey, Brainard, Lee,
        % et al., Vision Research, 2000, page 1808 bottom right.
        % (http://www.cns.nyu.edu/~tony/vns/readings/dacey-etal-2000.pdf)
        %
        % Also see Boycott & Wassle, 1991,  (European Journal of Neuroscience), Table 1.
        %
        % Likely the larger RF sizes measured physiological (Dacey et al.)
        % vs anatomically (B&W) reflect spread of signals among cones (via direct
        % gap junctions) and probably more important among cone bipolars
        % (via gap junctions with AII amacrine cells). - Fred
        
        %         @JRG incorporate preferntial cone selection to bipolar rfs
        switch obj.cellType
            case{'onDiffuse','offDiffuse','onDiffuseSBC'}
                % ecc = 0 mm yields 2x2 cone input to bp
                % ecc = 30 mm yields 5x5 cone input to bp
                sizeScale = floor(2 + 3/10*(p.Results.ecc)); 
                
                bpSizeCenter = sizeScale;
                bpSizeSurround = sizeScale;
%                 obj.sRFcenter = ones(sizeScale).*sizeScale^2;
%                 obj.sRFsurround = ones(sizeScale).*sizeScale^2;
                rfCenterBig = 1;%fspecial('gaussian',[bpSizeCenter,bpSizeCenter],1); % convolutional for now
                rfSurroundBig = 1;%fspecial('gaussian',[bpSizeSurround,bpSizeSurround],1); % convolutional for now
                
                obj.sRFcenter = rfCenterBig(:,:);
                obj.sRFsurround = rfSurroundBig(:,:);
                
            case{'onMidget','offMidget'}
                
                % ecc = 0 mm yields 1x1 cone input to bp
                % ecc = 30 mm yields 3x3 cone input to bp
                sizeScale = floor(1 + (2/10)*(p.Results.ecc)); 
                
                bpSizeCenter = sizeScale;
                bpSizeSurround = sizeScale;
                
%                 obj.sRFcenter = ones(sizeScale)./sizeScale^2;
%                 obj.sRFsurround = ones(sizeScale)./sizeScale^2;
                obj.sRFcenter = 1;%fspecial('gaussian',[bpSizeCenter,bpSizeCenter],1); % convolutional for now
                obj.sRFsurround = 1;%fspecial('gaussian',[bpSizeSurround,bpSizeSurround],1); % convolutional for now
             
        end
        if isfield(p.Results,'cellLocation')
            obj.cellLocation = p.Results.cellLocation;
        end
    end
    
    % set function, see bipolarSet for details
    function obj = set(obj, varargin)
        bipolarSet(obj, varargin{:});
    end
    
    % get function, see bipolarGet for details
    function val = get(obj, varargin)
        val = bipolarGet(obj, varargin{:});
    end
    
    % compute function, see bipolarCompute for details
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