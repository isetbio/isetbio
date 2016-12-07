classdef coneMosaic < hiddenHandle
    % Create a cone mosaic class
    %
    %   cMosaic =  coneMosaic('pigment', photoPigment, 'os', os);
    %
    % The cone mosaic defines the absorptions and photocurrent in an array
    % of cones. The default cone mosaic is rectangular.  There is a
    % subclass of coneMosaicHex (hexagonal sampling). That is implemented
    % by using a finer grid and placing the cones on a submosaic within the
    % larger rectangular mosaic.
    %
    % The cone quantum efficiencies are defined by both the macular pigment
    % and the pigment. The lens is incoporated as part of the optics.
    %
    % The isomerizations (absorptions, R*) are calculated using
    % coneMosaic.compute.
    %
    % The absorptions are converted to photocurrent by the method defined
    % in the outerSegment class, coneMosaic.os.  The compute for current is
    % coneMosaic.computeCurrent.
    %
    % HJ/JRG/BW ISETBIO Team, 2016
    
    % Keep the format for the comments.  This is used in doc coneMosaic
    properties (GetAccess=public, SetAccess=public) 

        name                % the name of the object
        
        pigment;            % Cone photopigment class
        macular;            % Macular class
        os;                 % outerSegment class

        absorptions;        % The spatial array of cone absorptions
                            % must be consistent with pattern
        current;            % The (x,y,t) of photocurrent, stored in os

        center;             % Center position of patch (x,y - meters)
        
        pattern;            % Pattern of K-LMS cones in the mosaic%
                            % Defines rows and cols, too 
        
        patternSampleSize;  % Separation between K-LMS pattern samples 
                            % For rectangular grid mosaics, this is set to
                            % the pigment width/height, i.e., the actual
                            % cone separation; For hexagonal grid mosaics
                            % (instances of the coneMosaicHex class), this
                            % is the separation between the rect grid nodes
                            % over which the lower resolution hex grid is
                            % sampled (NC)
        
        integrationTime;    % Cone temporal integration time (secs)
        emPositions;        % Eye movement positions (spatial units are number of cones).
                            % The number of positions controls number of
                            % frames to be computed
        noiseFlag;          % Absorption calculation noise (usually photon noise is on)
                            
        hdl;                % handle of the coneMosaic window
    end
    
    
    properties (Dependent)
        % Dependency shown in parenthesis
        
        wave;           % Wavelength samples (pigment)
        
        rows;           % number of rows in the cone mosaic (pattern)
        cols;           % number of cols in the cone mosaic (pattern)
        mosaicSize;     % [rows, cols] (pattern)
        patternSupport; % [X(:) Y(:)] axes (Pattern)
        
        width;          % width of cone mosaic in meters (patternSampleSize)
        height;         % height of cone mosaic in meters (patternSampleSize)
        fov;            % horizontal/vertical field of view assuming inf 
                        % scene distance and 17mm optics focal length
                        % (patternSampleSize, via height and width)
        tSamples        % Number of temporal samples
        
        coneLocs;       % cone locations in meters (pigment)
        qe;             % absorptance with macular pigment, but not lens
                        % which is accounted for in the oi representation.
                        % (pigment and macular)
        
        spatialDensity; % spatial density (ratio) of the K-LMS cones
        
        % absorptionsTimeAxis; % Time samples (absorptions and integrationTime) 
    end
    
    
    properties (Access=private)
        % spatial density (ratio) of the K-LMS cones
        %
        % There are two properties about this LMS ratio: sptailDensity_ and
        % spatialDensity. spatialDensity_ is a private variable of the
        % class while spatialDensity is a dependent variable.
        %
        % When setting the spatialDensity_ inside the class, we just change
        % this private property directly. When setting spatialDensity from
        % outside, the set.spatialDensity is executed and the cone mosaic
        % is regenerated.
        %
        % It's possible to have set function for spatialDensity directly.
        % However, the current Matlab does not allow update other
        % properties (e.g. mosaic pattern) in the set function for
        % non-dependent properties.
        %
        % HJ
        spatialDensity_;
    end
    
    properties (Constant)
        validNoiseFlags = {'random','none','frozen'};
    end
    
    methods    
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %   cMosaic =  coneMosaic('cone',cone,'os','os);
            
            % TODO:  We should take eccentricity and angle as an input
            % parameter and create cones of the appropriate size for that
            % retinal location
            
            % parse input
            p = inputParser;
            
            p.addParameter('name', 'cone mosaic', @ischar);
            
            % Included objects
            p.addParameter('pigment', photoPigment(), ...
                @(x) isa(x, 'photoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('os', osLinear(), @(x)(isa(x,'outerSegment')));

            % Mosaic parameters
            p.addParameter('center',[0 0], @(x)(numel(x) ==2));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('integrationTime', 0.05, @isscalar);
            
            % Computational features
            p.addParameter('emPositions', [0 0], @isnumeric);
            
            % How we handle coneMosaic noise
            vFunc = @(x)(ismember(lower(x), coneMosaic.validNoiseFlags));
            p.addParameter('noiseFlag', 'random', vFunc);
            
            p.parse(varargin{:});
            
            % set properties
            obj.name    = p.Results.name;
            obj.pigment = p.Results.pigment;
            obj.macular = p.Results.macular;
            obj.os      = p.Results.os;
            
            obj.center = p.Results.center(:)';
            obj.wave   = p.Results.wave;
            obj.spatialDensity_ = p.Results.spatialDensity(:);
            obj.integrationTime = p.Results.integrationTime;
            
            obj.noiseFlag = p.Results.noiseFlag;
            obj.emPositions = p.Results.emPositions;
            
            % Set the cone spacing and aperture given its eccentricity and
            % angle.  We could specify eye, but are we really sure about
            % the left right thing in human?
            % 
            % Units of returns are meters
            ecc = sqrt(sum(obj.center.^2));
            ang = atan2d(obj.center(2),obj.center(1));
            [spacing, aperture] = coneSize(ecc,ang);

            obj.pigment.pdWidth  = aperture;
            obj.pigment.pdHeight = aperture;
            obj.pigment.height = spacing;
            obj.pigment.width  = spacing;
            
            % See description of this parameter on the wiki page at
            % https://github.com/isetbio/isetbio/wiki/Cone-mosaic#note-on-the-patternsamplesize-variable
            obj.patternSampleSize = [obj.pigment.width, obj.pigment.height];
            
            % generate human cone mosaic pattern if not specified
            if isempty(p.Results.pattern)
                [~, obj.pattern] = humanConeMosaic(p.Results.size, ...
                    obj.spatialDensity_, obj.patternSampleSize(1));
            else
                obj.pattern = p.Results.pattern;
            end
            
            % Initialize the mosaic properties
            % obj.os.timeStep = obj.sampleTime;
            obj.os.patchSize = obj.width;
            
            % initialize listener
            % these listeners make sure the wavelength samples
            % in obj.pigment and obj.macular match
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);

        end
        
        
        %% get methods for dependent variables
        % http://www.mathworks.com/help/matlab/matlab_oop/specifying-methods-and-functions.html#bu4wzba
        % All functions that use dots in their names must be defined in the
        % classdef file, including:
        %   -  Converter methods that must use the package name as part of
        %      the class name because the class is contained in packages 
        %   -  Property set and get access methods
        function val = get.wave(obj)
            val = obj.pigment.wave;
        end
        
        function val = get.rows(obj)  % number of rows
            val = size(obj.pattern, 1);
        end
        
        function val = get.cols(obj)  % number of cols
            val = size(obj.pattern, 2);
        end
        
        function val = get.mosaicSize(obj)
            val = [obj.rows obj.cols];
        end
        
        function val = get.patternSupport(obj)
            x = (1:obj.cols) * obj.patternSampleSize(1); x = x - mean(x);
            y = (1:obj.rows) * obj.patternSampleSize(2); y = y - mean(y);
            [xx, yy] = meshgrid(x,y);
            val(:,:,1) = xx;
            val(:,:,2) = yy;
        end
        
        function val = get.width(obj)  % width of cone mosaic in meters
            val = obj.cols * obj.patternSampleSize(1);
        end
        
        function val = get.height(obj)  % height of cone mosaic in meters
            val = obj.rows * obj.patternSampleSize(2);
        end
        
        function val = get.fov(obj)  % horizontal/vertical field of view
            val = 2 * atand([obj.width obj.height]/2/0.017);
        end
        
        function val = get.coneLocs(obj) % cone locations in meters
            %
            x = (1:obj.cols) * obj.pigment.width; x = x - mean(x) + obj.center(1);
            y = (1:obj.rows) * obj.pigment.height; y = y - mean(y) + obj.center(2);
            
            [X, Y] = meshgrid(x, y);
            val = [X(:) Y(:)];
        end
        
        function val = get.qe(obj)
            % cMosaic.qe
            % compute effective absorptance with macular pigments
            %
            % The cone mosaic quantum efficiency is the product of the cone
            % photopigment absorptance times the macular pigment
            % transmittance times the cone photopigment peak absorption.
            %
            % The eye quantum efficiency, called in plot, includes the lens
            % transmittance.
            val = bsxfun(@times, obj.pigment.absorptance, ...
               obj.macular.transmittance)*diag(obj.pigment.peakEfficiency);
        end
        
        function val = get.spatialDensity(obj)
            val = obj.spatialDensity_;
        end
        
        function val = get.absorptions(obj)           
            val = double(obj.absorptions);
        end
        
        function val = get.current(obj)
            % Shouldn't have to do this any more, right?
            if isempty(obj.current), val = [];
            else                     val = double(obj.current);
            end
        end
        
        function val = get.tSamples(obj)
            % Computed from the 1st dimension of eye movement positions
            if isempty(obj.emPositions), val = 1; return; end
            
            % The number of eye positions defines the number of samples
            val = size(obj.emPositions,1);

            % If you want the number of time steps in sec as an axis:
            %
            %   val = (0:1:(tSamps-1)) * obj.integrationTime;

        end        
        
        %% set method for class properties
        function set.spatialDensity(obj, val)
            if all(obj.spatialDensity_(:) == val(:)), return; end
            obj.spatialDensity_ = val;
            [~, obj.pattern] = humanConeMosaic(obj.mosaicSize, ...
                    val, obj.patternSampleSize(1));
            obj.clearData();
        end
        
        function set.wave(obj, val)
            obj.pigment.wave = val(:);
            obj.macular.wave = val(:);
        end
       
        function set.absorptions(obj, val)
            obj.absorptions = single(val);
        end
        
        function set.current(obj, val)
            obj.current = single(val);
        end
        
        function set.fov(obj, val) % set field of view
            obj.setSizeToFOV(val);
        end
        
        function set.mosaicSize(obj, val)
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.spatialDensity_, obj.patternSampleSize(1));
                obj.os.patchSize = obj.width;
                obj.clearData();
            end
        end
        
        function set.rows(obj, val)
            obj.mosaicSize = [val obj.cols];
        end
        
        function set.cols(obj, val)
            obj.mosaicSize = [obj.rows val];
        end
        
        function set.noiseFlag(obj, val)
            if ischar(val) && (ismember(lower(val), coneMosaic.validNoiseFlags))
                obj.noiseFlag = val;
            else
                validNoiseFlags = coneMosaic.validNoiseFlags
                illegalValueForNoiseFlag = val
                error('This is an invalid value for coneMosaic.noiseFlag');
            end
        end
    end
    
    methods (Access=public)
        % Declare the compute method
        [absorptions, current, varargout] = compute(obj, oi, varargin);
        
        % Declare the compute method for a sequence of optical images viewed sequentially
        [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence(obj, oiSequence,  varargin)
        
        % Method returning the demosaiced isomerization maps and the corresponding sRGB rendition
        function [demosaicedAbsorptionsMap, sRGB] = demosaicedIsomerizationMaps(obj, varargin)
            [demosaicedAbsorptionsMap, sRGB] = obj.demosaicedResponses();
        end
        
        % Method returning the demosaiced photocurrent maps
        function demosaicedCurrentsMap = demosaicedPhotoCurrentMaps(obj, current)
            demosaicedCurrentsMap = obj.demosaicedResponses(current);
        end
        
        % Demosaicing method
        varargout = demosaicedResponses(obj, varargin);

        % Method that low-passes a mosaic's response using different space
        % constants for each of the L,M, and S-cone submosaic.
        [lowPassedResponse,  Lmap, Mmap, Smap] = lowPassActivationMap(response, theMosaic, spaceConstants);
        
        % Declare the computeCurrent method
        interpFilters = computeCurrent(obj, varargin);

        function val = timeAxis(obj)
            % Call:  obj.timeAxis;
            val = (0:1:(obj.tSamples-1)) * obj.integrationTime;
        end
        
    end

    methods (Static)
        [noisyImage, theNoise] = photonNoise(absorptions,varargin);
        resampledAbsorptionsSequence = tResample(absorptionsSequence, pattern, originalTimeAxis, resampledTimeAxis);
    end

    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        cpObj = copyElement(obj);
    end
    
    methods (Access = public, Hidden)
        absorptions = computeSingleFrame(obj, oi, varargin);
        absorptions = applyEMPath(obj, LMS, varargin);
    end        
    
    methods (Access = private)
        setWave(obj, src, ~)
    end

end