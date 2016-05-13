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
    patchSize;                       % size of retinal patch from sensor
    timeStep;                        % time step of simulation from sensor
    sRFcenter;                       % spatial RF of the center on the receptor grid
    sRFsurround;                     % spatial RF of the surround on the receptor grid
    temporalDifferentiator;          % differentiator function
    responseCenter;                  % Store the linear response of the center after convolution
    responseSurround;                % Store the linear response of the surround after convolution

end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

% Public methods
methods
    
    % Constructor
    function obj = bipolar(os)
        
        obj.patchSize = osGet(os,'patchSize');
        obj.timeStep = osGet(os,'timeStep');
        
        % Build spatial receptive field
        obj.sRFcenter = fspecial('gaussian',[2,2],1); % convolutional for now
        obj.sRFsurround = fspecial('gaussian',[2,2],1); % convolutional for now
        
        % Weight are caclulated offline by optimizing for the minimum error
        % between the differentiator signal and the bipolar output
        % calculated with filter from the deconvolution operation.        
        switch class(os)
            case{'osLinear'}
                coneW = -0.048515736811122e4; coneDiffW = -3.912753277547210e4;
            otherwise
                coneW = -0.1373; coneDiffW = -9.5355;
        end
        obj.temporalDifferentiator = @(x) coneW*x(:,2:end) + coneDiffW*diff(x,1,2);
        
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
