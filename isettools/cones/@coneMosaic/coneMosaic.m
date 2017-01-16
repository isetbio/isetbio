classdef coneMosaic < hiddenHandle
    %CONEMOSAIC - Create a cone mosaic object
    %   An object of the cone mosaic class defines parameters of the
    %   retinal cone mosaic, and enables computation of cone isomerizations
    %   and photocurrent response in an array of cones.
    %
    %   cMosaic =  coneMosaic(...,'PARAM1',val1,'PARAM2',val2,...) creates
    %   the cone mosaic object. Optional parameter name/value pairs are
    %   listed below.
    %
    %   The cone mosaic defines the absorptions and photocurrent in an array
    %   of cones. The default cone mosaic is rectangular.  There is a
    %   subclass CONEMOSAICHEX that allows more realistic hexagonal and
    %   variable spacing mosaics. That subclass is implemented by using a
    %   finer grid and placing the cones at a subset of locations on the
    %   grid.
    %
    %   The cone quantum efficiencies are defined by both the macular pigment
    %   and the pigment. The lens is incoporated as part of the optics.
    %
    %   The isomerizations (absorptions, R*) are calculated using
    %   coneMosaic.compute.
    %
    %   The isomerizations are converted to photocurrent by the method defined
    %   in the outerSegment class, coneMosaic.os.  The compute for current is
    %   coneMosaic.computeCurrent.
    %
    %   Optional parameter name/value pairs chosen from the following:
    %
    %   'name'            Mosaic name
    %   'pigment'         Cone photopigment object (defaults to what photoPigment() sets up).
    %   'macular'         Macular pigment object (defaults to what Macular() sets up).
    %   'os'              Outer segment object (defauls to what osLinear() sets up).
    %   'center'          Position of center of mosaic on the retina. Vector (default [0,0]).
    %   'wave'            Wavelength samples in nm. Vector (default 400:10:700).
    %   'pattern'         Cone type at each position (1-4, K,L,M,S) Matrix (default []).
    %   'spatialDensity'  Relative density of cone types, K,L,M,S. Vector (default [0 0.6 0.3 0.1]).
    %   'size'            Spatial size of mosaic (number of rows/cols). Vector (default [72 88]).
    %   'integrationTime' Value (default 0.005). Temporal integration in
    %                     sec. Keep this under 25 ms (0.025) if you are
    %                     computing photocurrent for decent numerical
    %                     accuracy.
    %   'emPositions'     Eye movement positions. Nx2 matrix (default [0 0] is
    %                     middle or cone mosaic, 1 unit is 1 cone for rect
    %                     (HEX?)
    %   'apertureBlur'    Blur by cone aperture? true/false (default false).
    %   'noiseFlag'       Add photon noise (default) or not. String (default
    %                     'random').  Valid values are 'random', 'frozen', or 'none'.
    %
    %  ISETBIO wiki: <a href="matlab:
    %  web('https://github.com/isetbio/isetbio/wiki/Cone-mosaic','-browser')">cone mosaic</a>.
    %
    % HJ/JRG/BW ISETBIO Team, 2016
    %
    % See also CONEMOSAICHEX, PHOTOPIGMENT, MACULAR, LENS, OUTERSEGMENT

    % Keep the format of the comments here.  They are used by the doc
    % command to produce useful documentation of the properites.
    properties (GetAccess=public, SetAccess=public)
        %NAME  The name of the object
        name;
        
        %PIGMENT  Cone photopigment object for mosaic 
        pigment;
        
        %MACULAR  Macular pigment object for mosaic
        macular;
        
        %OS  Outer segment object for mosaic
        os;
        
        %ABSORPTIONS  The spatial array of cone absorptions 
        %    This must be kept consistent with the mosaic pattern.
        absorptions;
        
        %CONEDARKNOISERATE  Mean rate of dark (thermal) isomerizations for each cone class.
        %    Expressed as a three-dimensional row vector,
        %    isomerizations/sec.
        %    Default values are [0 0 0] for backwards compatibility.
        %    Reasonable values are about [300 300 300].
        coneDarkNoiseRate;
        
        %CURRENT  The (x,y,t) of photocurrent 
        %    There is a comment that this is actually stored in the OS.
        %    Should check and explain this more.
        current;    
        
        %CENTER  Center position of patch (x,y - meters)
        %    This tells us where the mosaic is relative to the optical
        %    image.
        center;
        
        %PATTERN Pattern of KLMS cones in the mosaic
        %    K = 1, L = 2, M = 3, S = 4 (K means blank,no cone at that position)
        %    This defines the row and column dimensions of the grid on
        %    which the mosic is defined. 
        pattern;
        
        %PATTERNSAMPLESIZE  Separation between KLMS pattern samples
        %    For rectangular grid mosaics, this is set to the width/heigh field
        %    of the PIGMENT object, i.e., the actual cone separation.
        %
        %    For hexagonal grid mosaics (instances of the coneMosaicHex class), this is
        %    the separation between the rect grid nodes on which the
        %    cone positions are sampled.
        patternSampleSize; 
        
        %INTEGRATIONTIME  Cone temporal integration time (secs).
        %    Keep this under 25 ms (0.025) if you are computing
        %    photocurrent and want reasonable numerical accuracy.
        integrationTime;
        
        %EMPOSITIONS  Eye movement positions
        %    Spatial units are those of rectangular grid on which cones are
        %    specified.
        %
        %    The number of positions controls number of frames to be
        %    computed
        emPositions;
        
        %NOISEFLAG  Add noise to isomerizations?
        noiseFlag;
        
        %APERTUREBLUR  Blur by cone aperture when computing isomerizations (true/false)?
        apertureBlur;
        
        %HDL  Handle of the CONEMOSAIC window
        %    
        hdl
    end
    
    properties (Dependent)        
        %WAVE  Wavelength samples
        %    Depends on wavelength sampling in the PHOTOPIGMENT object
        %    specified in the PIGMENT property.
        wave;
        
        %ROWS  Number of rows in the cone mosaic
        %    Depends on size of PATTERN property.
        rows;  
        
        %COLS   Number of cols in the cone mosaic
        %    Depends on size of PATTERN property.
        cols;
        
        %MOSAICSIZE  Vector containing [rows cols]
        %    Depends on size of PATTERN property via ROWS and COLS.
        mosaicSize;
        
        %PATTERNSUPPORT  Matrix giving x,y positions of the underlying grid.
        %    In form [x(:) y(:)] in units of meters.  These are positions
        %    on the retina, relative to the center of the specified mosaic.
        %    It does not take the CENTER property into account.  It is
        %    based on the PATTERN and PATTERNSAMPLESIZE properties.
        patternSupport;
        
        %WIDTH  Width of cone mosaic in meters
        %    Depends on property patternSampleSize.
        width; 
        
        %HEIGHT  Height of cone mosaic in meters
        %    Depends on property patternSampleSize.
        height;
        
        %FOV  Vector containing nominal [hfov vfov] (field of view) in degrees  
        %    Computed assuming infiite scene distance and 17mm optics focal length
        %
        %    Depends on patternSampleSize via height and width properties.
        fov;
        
        %TSAMPLES  Number of temporal samples  
        tSamples
        
        %CONELOCS  Matrix giving x,y positions of the cone locations on the grid
        %    In form [x(:) y(:)] in units of meters.  These are positions
        %    on the retina, relative to the center of the specified mosaic,
        %    and these are offset by the CENTER property. It is based on
        %    the PATTERN and PIGMENT properties, as it uses the height and
        %    width specfied in the pigment rather than the mosaic to
        %    compute positions.
        coneLocs;
        
        %QE  Cone quantal efficiency (absorptance)
        %    The cone mosaic quantum efficiency is the product of the cone
        %    photopigment absorptance times the macular pigment
        %    transmittance times the cone photopigment peak efficientcy.
        %
        %    The quantum efficiency referred all the way to the cornea and
        %    used in the plot function, also includes the lens
        %    transmittance. The lens transmittance is specified by the LENS
        %    object in the oi representation.
        %
        %    The QE depends on the PIGMENT and MACULAR properties.
        qe;
        
        %SPATIALDENSITY  Spatial density of the KLMS cones  
        spatialDensity; 
    end
    
    % properties (Constant)
    %     % dictionary with all coneMosaic-specific warning labels 
    %     % and their default states (either 'on', or 'off')
    %     warnings = containers.Map(...
    %         {...
    %             'ISETBIO:ConeMosaic:computeForOISequence:displaySizeInfo', ...    % whether computeForOISequence() should display the size of the absorptions matrix
    %             'ISETBIO:ConeMosaic:osCompute:displayFrozenNoiseSeed' ...         % whether osCompute() should display the seed for the frozen noise
    %             'ISETBIO:ConeMosaic:integrationTimeSetToGreaterThan25msec' ...    % whether setting coneMosaic.integrationTime > 25 ms should display a warning
    %         }, ...
    %         { ...
    %             'off', ...                % default state for 'ISETBIO:ConeMosaic:computeForOISequence:displaySizeInfo': 
    %             'off' ...                 % default state for 'ISETBIO:ConeMosaic:osCompute:displayFrozenNoiseSeed'
    %             'off' ...                 % default state for 'ISETBIO:ConeMosaic:integrationTimeSetToGreaterThan25msec'
    %         });
    % end
    
    properties (Access=private)
        %SPATIALDENSITY_  Ratio of KLMS cones used to generate pattern
        %
        %    There are two properties about this ratio: sptiallDensity_ and
        %    spatialDensity. spatialDensity_ is a private variable of the
        %    class while spatialDensity is a dependent variable.
        %
        %    When setting the spatialDensity_ inside the class, we just change
        %    this private property directly. When setting spatialDensity from
        %    outside, the set.spatialDensity is executed and the cone mosaic
        %    is regenerated.
        %
        %    It's possible to have set function for spatialDensity directly.
        %    However, the current Matlab does not allow update other
        %    properties (e.g. mosaic pattern) in the set function for
        %    non-dependent properties.
        spatialDensity_;
    end
    
    properties (Constant)
        %VALIDNOISEFLAGS  Cell array of strings containing valid values for noise flags.
        validNoiseFlags = {'none','frozen','random'};
    end
    
    methods
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %   cMosaic =  coneMosaic('cone',cone,'os','os);
            %
            % TODO:  We should take eccentricity and angle as an input
            % parameter and create cones of the appropriate size for that
            % retinal location
            p = inputParser;
            p.addParameter('name', 'cone mosaic', @ischar);
            p.addParameter('pigment', photoPigment(),@(x) isa(x, 'photoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));  
            p.addParameter('os', osLinear(), @(x)(isa(x,'outerSegment'))); 
            p.addParameter('center',[0 0], @(x)(numel(x) ==2));   
            p.addParameter('wave', 400:10:700, @isnumeric);       
            p.addParameter('pattern', [], @isnumeric);            
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);  
            p.addParameter('size', [72 88], @isnumeric);         
            p.addParameter('integrationTime', 0.005, @isscalar);  
            p.addParameter('emPositions', [0 0], @isnumeric);
            p.addParameter('apertureBlur', false, @islogical);
            p.addParameter('noiseFlag', 'random', @(x)(ismember(lower(x), coneMosaic.validNoiseFlags)));
            p.parse(varargin{:});
            
            % Set properties
            obj.name    = p.Results.name;
            obj.pigment = p.Results.pigment;
            obj.macular = p.Results.macular;
            obj.os      = p.Results.os;
            
            obj.center = p.Results.center(:)';
            obj.wave   = p.Results.wave;
            obj.spatialDensity_ = p.Results.spatialDensity(:);
            obj.integrationTime = p.Results.integrationTime;
            
            obj.coneDarkNoiseRate = [0 0 0];
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
            
            % Blur by aperture
            obj.apertureBlur = p.Results.apertureBlur;
            
            % Initialize listener
            %
            % These listeners make sure the wavelength samples
            % in obj.pigment and obj.macular match
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);     
            
            % Set default warnings
            % coneMosaic.allWarnings('default');
        end
        
        %% Get methods for dependent variables
        % http://www.mathworks.com/help/matlab/matlab_oop/specifying-methods-and-functions.html#bu4wzba
        % All functions that use dots in their names must be defined in the
        % classdef file, including:
        %   -  Converter methods that must use the package name as part of
        %      the class name because the class is contained in packages
        %   -  Property set and get access methods
        
        function val = get.wave(obj)
            val = obj.pigment.wave;
        end
        
        function val = get.rows(obj)
            val = size(obj.pattern, 1);
        end
        
        function val = get.cols(obj)
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
        
        function val = get.width(obj)
            val = obj.cols * obj.patternSampleSize(1);
        end
        
        function val = get.height(obj)
            val = obj.rows * obj.patternSampleSize(2);
        end
        
        function val = get.fov(obj)
            val = 2 * atand([obj.width obj.height]/2/0.017);
        end
        
        function val = get.coneLocs(obj) % cone locations in meters
            x = (1:obj.cols) * obj.pigment.width; x = x - mean(x) + obj.center(1);
            y = (1:obj.rows) * obj.pigment.height; y = y - mean(y) + obj.center(2);
            
            [X, Y] = meshgrid(x, y);
            val = [X(:) Y(:)];
        end
        
        function val = get.qe(obj)
            % obj.qe  Compute effective absorptance including macular pigment
            val = bsxfun(@times, obj.pigment.absorptance, ...
                obj.macular.transmittance)*diag(obj.pigment.peakEfficiency);
        end
        
        function val = get.spatialDensity(obj)
            val = obj.spatialDensity_;
        end
        
        % Shouldn't have to do this any more, right?
        function val = get.current(obj)
            if isempty(obj.current)
                val = [];
            else
                val = double(obj.current);
            end
        end
        
        function val = get.absorptions(obj)
            val = double(obj.absorptions);
        end
        
        function val = get.tSamples(obj)
            % Computed from the 1st dimension of eye movement positions
            if isempty(obj.emPositions), val = 1; return; end
            
            % The number of eye positions defines the number of samples
            val = size(obj.emPositions,1);
        end
        
        %% Set methods for class properties
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
        
        function set.integrationTime(obj, val)
            % if (val > 25/1000)
            %     warning('ISETBIO:ConeMosaic:integrationTimeSetToGreaterThan25msec', ...
            %         'Setting the coneMosaic.integrationTime > 25ms is not recommended if you are interested in photocurrent computations. Assigned value: %2.2fmsec.', val*1000);
            % end
            obj.integrationTime = val;
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
                s = sprintf('%s ', coneMosaic.validNoiseFlags{:});
                error('''%s'' is an invalid value for coneMosaic.noiseFlag. Choose one from: %s ', val,s);
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
        [interpFilters, meanCur] = computeCurrent(obj, varargin);
        
        function val = timeAxis(obj)
            % Call:  obj.timeAxis;
            val = (0:1:(obj.tSamples-1)) * obj.integrationTime;
        end
        
    end
    
    methods (Static)
        [noisyImage, theNoise] = photonNoise(absorptions,varargin);
        allWarnings(state);        % valid state values : 'default', 'on', 'off'
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