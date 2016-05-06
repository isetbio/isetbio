classdef bipolar < handle
% Define the bipolar cell class.
% 
% The bipolar class allows the simulation of nonlinear subunits within
% retinal ganglion cell spatial receptive fields. 
% 
% Input: the cone response over time from an outer segment object.
% 
% Output: the bipolar response over time.
% 
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    
    % The type of computational model for the RGC spikes
    % model;
    
    cellLocation;
    patchSize;
    timeStep;
    sRF;           % spatial RF of the center on the receptor grid
    tIR;            % temporal impulse response of the center
    threshold;     % threshold for nonlinear output
    response;      % Store the linear response after convolution
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
        obj.sRF = fspecial('gaussian',[7,7],1); % convolutional for now
        
        % Build temporal impulse response from M cone filter
        % bipolar temopral response = M cone filter + 0.5 * d/dt (M cone filter)
        osFilter = os.mConeFilter;        
        osFilterDerivativeShort = diff(osFilter);
        % osFilterDerivative = interp1(1:length(osFilterDerivativeShort),osFilterDerivativeShort,1:length(osFilter))';
        osFilterDerivative = [osFilterDerivativeShort; osFilterDerivativeShort(end)];
        % obj.tIR = osFilter + .5*osFilterDerivative;
        obj.tIR = 1;
        % Set threshold
        % obj.threshold = 0;
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
