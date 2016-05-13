classdef bipolar < handle
% Define the bipolar cell class
% 
% The bipolar class simulates processing from the cone outer segment
% current to the voltage output of the bipolar.
%
% Bipolar computation starts with the outer segment photocurrent, applies a
% linear temporal impulse response, and then weights over space.
% 
% The bipolar temporal impulse response (IR) is modeled using a widely used
% approximation from the literature (WHO?)
%
%   irBipolar(t) = w1*irCone(t) + w2 * (d/dt (irCone))
%
% where w1 is close to zero and w2 is a large negative number.  These
% values are stored in the bipolar object. We have default values for
% biophys and otherwise.
%
% The bipolar object sums over space using nonlinear subunits. 
% 
% Input: 
%   Cone outer segment object.  If none is passed in, we use the
%   osCreate('linear') object by default.
% 
% Output:
%  A bipolar object
%
% Example:
%   bp = bipolar;
%   bp = bipolar(osCreate('linear'));
%   bp = bipolar(osCreate('biophys'));
% 
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    
    % We are going to work on cellLocation.  At the moment, there is a
    % bipolar centered at every third cone in the os mosaic.  This is
    % implemented in the compute function.  We plan on specifying the
    % locations here, though.
    cellLocation;                    % location of bipolar RF centers
    sRFcenter;                       % spatial RF of the center on the receptor grid
    sRFsurround;                     % spatial RF of the surround on the receptor grid
    
    patchSize;                       % size of cone mosaic (os)
    timeStep;                        % time step of simulation from os
  
    responseCenter;                  % Store the linear response of the center after convolution
    responseSurround;                % Store the linear response of the surround after convolution

    temporalDifferentiator;          % differentiator function
    coneW;                           % Default weights for linear os case
    coneDiffW;
    
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

% Public methods
methods
    
    % Constructor
    function obj = bipolar(os)
        
        if ~exist('os','var'), os = osCreate('linear'); end
        
        obj.patchSize = osGet(os,'patchSize');
        obj.timeStep = osGet(os,'timeStep');
        
        % Build spatial receptive field
        obj.sRFcenter = fspecial('gaussian',[2,2],1); % convolutional for now
        obj.sRFsurround = fspecial('gaussian',[2,2],1); % convolutional for now
        
        switch class(os)
            case 'osBiophys'
                obj.coneW = -0.1373; obj.coneDiffW = -9.5355;
            otherwise
                % Use linear values
                obj.coneW     = -0.0485e4; 
                obj.coneDiffW = -3.912e4;
        end
        
        % Weight are caclulated offline by optimizing for the minimum error
        % between the differentiator signal and the bipolar output
        % calculated with filter from the deconvolution operation.        
        obj.temporalDifferentiator = @(x) obj.coneW*x(:,2:end) + obj.coneDiffW*diff(x,1,2);
        
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
        % Requires an outer segment with a time series
        % Returns the voltage over time for each cell in the bipolar mosaic
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
