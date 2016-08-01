classdef ir < handle
    %% The inner retina class stores general properties of the retinal patch 
    % and stores the rgcMosaic objects in its mosaic property field.
    %
    %   obj = ir(inputObj, params);    [usually called internally from irCreate]
    %
    % An ir object takes as input a bipolar object or an outerSegment object.
    % The ir (inner retina) object stores basic properties about the inner
    % retina such as the position of the simulated retinal patch.
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
    % Methods: set, get, compute, plot
    %       see individual m-files for details.
    %
    % Examples:
    %
    %   innerRetina1 = irCreate(bp,'name','myRGC');
    %
    %   params.name = 'Macaque inner retina 1'; %
    %   innerRetina2 = ir(bp, params);
    %
    % (c) isetbio team
    %
    % 9/2015 JRG
    % 7/2016 JRG updated
    
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
        function obj = ir(inputObj, varargin)
            
            % parse input
            p = inputParser;
            p.addRequired('inputObj');
            p.addParameter('eyeSide','left',@ischar);
            p.addParameter('eyeRadius',0,@isnumeric);
            p.addParameter('eyeAngle',0,@isnumeric);
            p.addParameter('name','ir1',@ischar);
            p.addParameter('species','macaque',@ischar);
            
            p.parse(inputObj,varargin{:});
            
            obj.eyeSide   = p.Results.eyeSide;
            obj.eyeRadius = p.Results.eyeRadius;
            obj.eyeAngle  = p.Results.eyeAngle;
            obj.name      = p.Results.name;            
            
            switch class(inputObj)
                case{'osDisplayRGB'}
                    obj.spacing = osGet(inputObj,'patch size'); % Cone width
                    obj.timing  = osGet(inputObj,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(osGet(inputObj,'rgbData'));
                case{'osIdentity'}
                    obj.spacing = osGet(inputObj,'patch size'); % Cone width
                    obj.timing  = osGet(inputObj,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(osGet(inputObj,'photonRate'));
                case{'bipolar'}
                    obj.spacing = bipolarGet(inputObj,'patch size'); % Bipolar width
                    obj.timing  = bipolarGet(inputObj,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~, ~] = size(bipolarGet(inputObj,'bipolarResponseCenter'));
                otherwise
                    obj.spacing = osGet(inputObj,'patch size'); % Cone width
                    obj.timing  = osGet(inputObj,'time step'); % Temporal sampling
                    [obj.row, obj.col, ~] = size(osGet(inputObj,'coneCurrentSignal'));
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
        function obj = compute(obj, inputObj, varargin)
            obj = irCompute(obj,  inputObj, varargin{:});
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
        
        % normalize function, see irNormalize
        function obj = normalize(obj,varargin)
            obj = irNormalize(obj, varargin{:});
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
    end
    
end


