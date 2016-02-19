classdef ir < handle
%% The inner retina class represents the mosaic of cells in the inner retina 
%
%   obj = ir(os, params);    %[usually called internally from irCreate]
% 
% An ir object is created from the isetbio @outerSegment object. The ir
% (inner retina) object stores basic properties about the species and the
% position of the simulated retinal patch. 
% 
% See Pillow, Jonathan W., et al. "Spatio-temporal correlations and visual
% signalling in a complete neuronal population." Nature 454.7207 (2008)
% and Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries 
% in ON and OFF ganglion cells of primate retina." The Journal of 
% Neuroscience 22.7 (2002).
% 
% Properties:
%    name:  animal, ir; example: 'macaque ir'
%    input: RGB stim or outersegment
%    temporalEquivEcc: calculated from retinal position, see retinalLocationToTEE        
%    mosaic: a cell array, where each cell is an rgcMosaic object, 
%               which is a subclass of the ir object.
% 
% Methods: intialize, set, get, compute, plot, movie
%       see individual .m files for details.
% 
% 
% Examples: 
% 
%   os  = osCreate('identity');
%   innerRetina1 = irCreate(os,'GLM','name','myRGC'); 
% 
%   params.name = 'Macaque inner retina 1'; % 
%   innerRetina2 = ir(os, params);
% 
% (c) isetbio
% 
% 9/2015 JRG
% 

%%
% Public read/write properties
properties
end

% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
    
end

% Protected properties; Methods of the parent class and all of its
% subclasses can set these.
properties (SetAccess = protected)
    % The type of computation is specified by the ir
    % subclass
    % These are the general inner retina parameters
    name;                  % Name of innerRetina, typically species name
    
    % Input spacing and count
    row;                   % N Stimulus row samples
    col;                   % N Stimulus col samples
    spacing;               % Stimulus input spacing (um)
    timing;                % Stimulus temporal sampling (sec)
    
    % Eye properties
    eyeSide;               % Left or right eye
    eyeRadius;             % Position of patch in radius
    eyeAngle;              % and angle (degrees)
    temporalEquivEcc;      % Temporal equivalent eccentricity
    
    % The ganglion cell types as a cell array
    mosaic;
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

% Public methods
methods
    function obj = ir(os, params)
        obj.initialize(os, params);        
    end
    
    function obj = mosaicCreate(varargin)
        obj = rgcMosaicCreate(varargin{:});
    end
    
    % set function, see irSet
    function obj = set(obj, varargin)
        obj = irSet(obj, varargin{:});
    end
    
    % get function, see irGet
    function val = get(obj, varargin)
        val = irGet(obj, varargin{:});
    end
    
    % IR Compute functions, that loop over the rgc mosaics
    function obj = compute(obj, outerSegment, varargin)
        % Calls computeContinuous and then computeSpikes
        obj = irCompute(obj,  outerSegment, varargin{:});
    end
    
    function obj = computeContinuous(obj,varargin)
        obj = irComputeContinuous(obj,varargin{:});
    end
    
    function obj = computeSpikes(obj, varargin)
        obj = irComputeSpikes(obj,  varargin{:});
    end
    
    % plot function, see irPlot
    function plot(obj, varargin)
        irPlot(obj, varargin{:});
    end
    
    % movie function, see irMovie
    function movie(obj, outersegment, varargin)
        irMovie(obj, outersegment, varargin{:});
    end
    
end

% Methods that must only be implemented in the subclasses.
% If a subclass does not implement each and every of these methods
% it cannot instantiate objects.

methods (Abstract, Access=public)
end

% Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
    
    spConvolve(obj);
    fullConvolve(obj);
    
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
    initialize(obj, sensor, outersegment, varargin);
end

end


    