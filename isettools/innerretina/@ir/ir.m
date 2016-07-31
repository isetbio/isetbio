classdef ir < handle
    %% The inner retina class represents cell mosaics
    %
    %   obj = ir(os, params);    %[usually called internally from irCreate]
    %
    % An ir object takes as input a bipolar object or an outerSegment object.
    % The ir (inner retina) object stores basic properties about the inner
    % retina such as the position of the simulated retinal patch.
    %
    % GLM model
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
    %       see individual m-files for details.
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
    properties (SetAccess = public, GetAccess = public)
        
        % We typically run a single trial
        numberTrials = 1;
        
        % The ganglion cell types as a cell array
        mosaic;
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
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        function obj = ir(os, varargin)
            
            % parse input
            p = inputParser;
            p.addRequired('os');
            p.addParameter('eyeSide','left',@ischar);
            p.addParameter('eyeRadius',0,@isnumeric);
            p.addParameter('eyeAngle',0,@isnumeric);
            p.addParameter('name','ir1',@ischar);
            p.addParameter('species','macaque',@ischar);
            
            p.parse(os,varargin{:});
            
            obj.eyeSide   = p.Results.eyeSide;
            obj.eyeRadius = p.Results.eyeRadius;
            obj.eyeAngle  = p.Results.eyeAngle;
            obj.name      = p.Results.name;            
            
            switch class(os)
                case{'osDisplayRGB'}
                    obj.spacing = osGet(os,'patch size'); % Cone width
                    obj.timing  = osGet(os,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(osGet(os,'rgbData'));
                case{'osIdentity'}
                    obj.spacing = osGet(os,'patch size'); % Cone width
                    obj.timing  = osGet(os,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(osGet(os,'photonRate'));
                case{'bipolar'}
                    obj.spacing = bipolarGet(os,'patch size'); % Bipolar width
                    obj.timing  = bipolarGet(os,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(bipolarGet(os,'bipolarResponseCenter'));
                otherwise
                    obj.spacing = osGet(os,'patch size'); % Cone width
                    obj.timing  = osGet(os,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~] = size(osGet(os,'coneCurrentSignal'));
            end
            
            % Initialize the mosaic property but do not generate any mosaics
            obj.mosaic = cell(1); % populated by mosaicCreate method
            
            % Get the TEE.
            obj.temporalEquivEcc = retinalLocationToTEE(obj.eyeAngle, obj.eyeRadius, obj.eyeSide);
            
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
        
        function obj = computeLinearSTSeparable(obj,varargin)
            obj = irComputeLinearSTSeparable(obj,varargin{:});
        end
        
        function obj = computeSpikes(obj, varargin)
            obj = irComputeSpikes(obj,  varargin{:});
        end
        
        % plot function, see irPlot
        function plot(obj, varargin)
            irPlot(obj, varargin{:});
        end
        
        function obj = normalize(obj,varargin)
            obj = irNormalize(obj, varargin{:});
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
        timeConvolve(obj);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj, sensor, outersegment, varargin);
    end
    
end


