classdef coneMosaic < hiddenHandle
    % Create a cone mosaic object - to be deprecated and replaced by cMosaic
    %
    
    % Syntax:
    %   cMosaic = coneMosaic;
    %   cMosaic = coneMosaic('pigment', pp);
    %
    % Description:
    %    An object of the cone mosaic class defines parameters of the
    %    retinal cone mosaic, and enables computation of cone
    %    isomerizations and photocurrent response in an array of cones.
    %
    %    cMosaic =  coneMosaic(..., 'PARAM1', val1, 'PARAM2', val2, ...)
    %    creates the cone mosaic object. Optional parameter name/value
    %    pairs are listed below.
    %
    %    The cone mosaic defines the absorptions and photocurrent in an
    %    array of cones. The default cone mosaic is rectangular.  There is
    %    a subclass coneMosaicHex that allows more realistic hexagonal and
    %    variable spacing mosaics. That subclass is implemented by using a
    %    finer grid and placing the cones at a subset of locations on the
    %    grid.
    %
    %    The cone quantum efficiencies are defined by both the macular
    %    pigment and the pigment. The lens is incoporated as part of the
    %    optics.
    %
    %    The isomerizations (absorptions, R*) are calculated using method
    %    compute.
    %
    %    The isomerizations are converted to photocurrent by the method
    %    defined in the outerSegment class, coneMosaic.os.  The compute
    %    method for current is computeCurrent.
    %
    %    [DHB NOTE: Need to update coordinate system for global isetbio
    %    consistency and document fact that key/value pairs for underlying
    %    read routines get passed on. This will take a bit of careful work
    %    and coordination.  Probably needs to be done on a branch and
    %    tested.]
    %
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    cMosaic           The created coneMosaic object.
    %
    % Optional key/value pairs:
    %    'name'             - String. Mosaic name. Default 'cone mosaic'.
    %    'pigment'          - Cone photopigment object. Default is set by
    %                         calling photoPigment().
    %    'macular'          - Macular pigment object. Default is set by
    %                         calling Macular().
    %    'os'               - Outer segment object. Default is set by
    %                         calling osLinear().
    %    'center'           - Vector. Default [0 0]. Position of center of
    %                         mosaic on the retina. In meters.
    %    'wave'             - Vector. Default 400:10:700. Wavelength
    %                         samples in nm.
    %    'pattern'          - Matrix. Default []. Cone type at each
    %                         position (1-4, K, L, M, S).
    %    'spatialDensity'   - Vector. Default [0 0.6 0.3 0.1]. Relative
    %                         density of cone types, K, L, M, S.
    %    'size'             - Vector. Default [72 88]. Spatial size of
    %                         mosaic (number of rows/cols).
    %    'integrationTime'  - Value. Default 0.005. Temporal integration in
    %                         sec. Keep this under 25 ms (0.025) if you are
    %                         computing photocurrent for decent numerical
    %                         accuracy.
    %    'micronsPerDegree' - Microns/degree. Default 300.
    %    'emPositions'      - Eye movement positions. N x 2 matrix (default
    %                         [0 0] is middle of cone mosaic, 1 unit is 1
    %                         cone for rect.
    %    'apertureBlur'     - Boolean. Default false. Blur by cone aperture
    %    'noiseFlag'        - String. Default 'random'. Add photon noise
    %                         (default) or not. Valid values are 'random', 
    %                         'frozen', or 'none'.
    %    'eccentricityunits' - 'm' or 'deg'
    %    'useParfor'
    %
    % Notes:
    %    * TODO: We should take eccentricity and angle as an input
    %      parameter and create cones of the appropriate size for that
    %      retinal location. (taken from constructor)
    %    * TODO: Allow units to be specified as a part of varargin for
    %      size() function below.
    %
    % References:
    %    ISETBIO wiki: https://github.com/isetbio/isetbio/wiki/Cone-mosaic
    %
    % See Also:
    %    coneMosaicHex, photoPigment, Macular, lens, outersegment
    %
    
    % History:
    %    xx/xx/16  HJ/JRG/BW  ISETBIO Team, 2016
    %    02/26/18  jnm        Formatting
    %    06/16/18  NPC        Function to determine whether a mosaic can
    %                         support ecc-based cone efficiency corrections
    %    12/18/20  dhb        Use new quantalEfficiency property of photoPigment.
    
    properties (GetAccess=public, SetAccess=public)
        %name - The name of the object
        name;

        % Species - {'human', 'macaque', 'mouse'}
        species;

        %pigment - Cone photopigment object for the mosaic
        pigment;

        %macular - Macular pigment object for mosaic
        macular;

        %os - Outer segment object for mosaic
        os;

        %absorptions - The spatial array of cone absorptions 
        %   Absorptions must be kept consistent with the mosaic pattern.
        absorptions;

        %coneDarkNoiseRate - Mean dark isomerizations rate per cone class.
        %   Dark is thermal
        %   Expressed as a three-dimensional row vector, 
        %   isomerizations/sec.
        %   Default values are [0 0 0] for backwards compatibility.
        %   Fred Rieke (personal communication with DHB) recommended [300 300 300].
        %   However, a lower value, e.g., [100 100 100] might be more
        %   appropriate according to Table 5 of 
        %   Koenig and Hofer (2011): The Absolute threshold of cone vision
        coneDarkNoiseRate;

        %current - The (x, y, t) of photocurrent
        %    There is a comment that this is actually stored in the OS.
        %    Should check and explain this more.
        current;

        %center - Center position of patch (x, y - meters)
        %   The mosaic center relative on the retina.  Used to set cone
        %   aperture and density parameters
        center;

        %whichEye - String ('left' or 'right') to indicate which eye
        whichEye;

        %pattern - Pattern of KLMS cones in the mosaic
        %   K = 1, L = 2, M = 3, S = 4
        %   (K means blank, no cone at that position)
        %   This defines the row and column dimensions of the grid on
        %   which the mosic is defined. 
        pattern;

        %patternSampleSize - Separation between KLMS pattern samples
        %   For rectangular grid mosaics, this is set to the width/height
        %   field of the PIGMENT object, i.e., the actual cone separation.
        %
        %   For hexagonal grid mosaics (instances of the coneMosaicHex
        %   class), this is the separation between the rect grid nodes on
        %   which the cone positions are sampled.
        patternSampleSize;

        %integrationTime - Cone temporal integration time (secs).
        %   Keep this under 25 ms (0.025) if you are computing
        %   photocurrent and want reasonable numerical accuracy.
        integrationTime;

        %micronsPerDegree - How many microns/degree. Default 300.
        %   Defaults to 300 which is appropriate for the central 
        %   region of the human eye.
        micronsPerDegree;

        %emPositions - Eye movement positions. Spatial and numerical.
        %   Spatial units are those of rectangular grid on which cones are
        %   specified. The number of positions controls number of frames to
        %   be computed
        emPositions;

        %noiseFlag - String. Add noise to isomerizations?
        noiseFlag;

        %apertureBlur - Boolean. Blur by cone ap. when computing iso's?
        %   The boolean indicating whether or not to blur by cone aperture
        %   when computing the isomerizations?
        apertureBlur;

        %useParfor - Boolean. Whether to use parfor to accelerate mosaic
        %    construction
        useParfor;
        
        %hdl  Handle of the CONEMOSAIC window
        hdl
    end

    properties (Dependent)
        %wave - Wavelength samples
        %   Depends on wavelength sampling in the PHOTOPIGMENT object
        %   specified in the PIGMENT property.
        wave;

        %rows - Number of rows in the cone mosaic
        %   Depends on size of PATTERN property.
        rows;

        %cols - Number of cols in the cone mosaic
        %   Depends on size of PATTERN property.
        cols;

        %mosaicSize - Vector containing [rows cols]
        %   Depends on size of PATTERN property via ROWS and COLS.
        mosaicSize;

        %patternSupport - Matrix giving x, y positions of underlying grid.
        %   In form [x(:) y(:)] in units of meters.  These are positions
        %   on the retina, relative to the center of the specified mosaic.
        %   It does not take the CENTER property into account.  It is
        %   based on the PATTERN and PATTERNSAMPLESIZE properties.
        patternSupport;

        %width - Width of cone mosaic in meters
        %   Depends on property patternSampleSize.
        width;

        %height - Height of cone mosaic in meters
        %   Depends on property patternSampleSize.
        height;

        %fov - Vector containing nominal [hfov vfov] FOV in degrees.
        %   The field of view is computed assuming infinite scene distance
        %   and 17mm optical focal length.
        %
        %   Depends on patternSampleSize via height and width properties.
        fov;

        %innerSegmentCoverage - The mosaic's inner segment coverage factor
        %   The retinal coverage factor of the mosaic computed as 
        %   (number of cones) x (inner segment area) / (mosaic area)
        %
        %   Depends on the inner segment area and the mosaic area.
        innerSegmentCoverage;
        
        %coverage - The  mosaic's coverage factor
        %   The retinal coverage factor of the mosaic computed as 
        %   (number of cones) x (cone area) / (mosaic area)
        %
        %   Depends on the cone area and the mosaic area.
        coverage;
        
        %tSamples  Number of temporal samples
        tSamples

        %coneLocs - Matrix giving x, y positions of the cone locs on grid
        %   In form [x(:) y(:)] in units of meters.  These are positions
        %   on the retina, relative to the center of the specified mosaic, 
        %   and these are offset by the CENTER property. It is based on
        %   the PATTERN and PIGMENT properties, as it uses the height and
        %   width specfied in the pigment rather than the mosaic to
        %   compute positions.
        coneLocs;

        %qe - Cone quantal efficiency (absorptance)
        %   The cone mosaic quantum efficiency is the product of the cone
        %   photopigment absorptance times the macular pigment
        %   transmittance times the cone photopigment peak efficientcy.
        %
        %   The quantum efficiency referred all the way to the cornea and
        %   used in the plot function, also includes the lens
        %   transmittance. The lens transmittance is specified by the LENS
        %   object in the oi representation.
        %
        %   The QE depends on the PIGMENT and MACULAR properties.
        qe;

        %spatialDensity - Spatial density of the KLMS cones
        spatialDensity;
    end

    properties (Access=private)
        %spatialDensity_ - Ratio of KLMS cones used to generate pattern
        %
        %   There are two properties about this ratio: sptiallDensity_ and
        %   spatialDensity. spatialDensity_ is a private variable of the
        %   class while spatialDensity is a dependent variable.
        %
        %   When setting the spatialDensity_ inside the class, we just
        %   change this private property directly. When setting
        %   spatialDensity from outside, the set.spatialDensity is executed
        %   and the cone mosaic is regenerated.
        %
        %   It's possible to have set function for spatialDensity directly.
        %   However, the current Matlab does not allow update other
        %   properties (e.g. mosaic pattern) in the set function for
        %   non-dependent properties.
        spatialDensity_;
    end
    
    properties (Constant)
        %validNoiseFlags - Cell string array containing the valid options.
        %   These options are: 'random', 'frozen', and 'none'.
        validNoiseFlags = {'none', 'frozen', 'random'};
    end

    methods
        % Constructor
        function [obj, cmParams] = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %
            % Syntax:
            %	[obj, cmParams] =  coneMosaic([varargin])
            %
            % Description:
            %    Initialize the cone mosaic class.
            %
            % Inputs:
            %    None required.
            %
            % Outputs:
            %    obj - The created cone mosaic object.
            %    cmParams - Struct with the params for this object
            %
            % Optional key/value pairs:
            %    Listed above.
            %
            % Notes:
            %    * TODO:  At present we take center (in meters) and compute
            %    eccentricity and angle. We use that to figure out cone
            %    aperture size and spacing. I wonder if we should take
            %    eccentricity and angle as an input parameter and computer
            %    center
            %
            %    * TODO:  Use ieParamFormat to simplify naming and allow
            %    spacing of the input arguments.
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('name', 'cone mosaic', @ischar);
            p.addParameter('pigment', photoPigment(), ...
                @(x) isa(x, 'photoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('os', [], ...
                @(x)(isempty(x) || isa(x, 'outerSegment')));
            p.addParameter('center', [0 0], @(x)(numel(x) == 2));
            
            % Should check for valid.  Maybe meters or degrees?  Microns?
            p.addParameter('eccentricityunits', 'm', @ischar); 

            p.addParameter('whichEye', 'left', ...
                @(x) ismember(x, {'left', 'right'}));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('micronsPerDegree', 300, @isnumeric);
            p.addParameter('integrationTime', 0.005, @isscalar);
            p.addParameter('emPositions', [0 0], @isnumeric);
            p.addParameter('apertureBlur', false, @islogical);
            p.addParameter('noiseFlag', 'random', ...
                @(x)(ismember(lower(x), coneMosaic.validNoiseFlags)));
            p.addParameter('useParfor', false, @islogical);
            p.parse(varargin{:});

            
            % Set properties
            obj.name    = p.Results.name;
            obj.pigment = p.Results.pigment;
            obj.macular = p.Results.macular;

            obj.center   = p.Results.center(:)';
            obj.whichEye = p.Results.whichEye;
            obj.wave     = p.Results.wave;
            obj.spatialDensity_   = p.Results.spatialDensity(:);
            obj.integrationTime   = p.Results.integrationTime;
            obj.micronsPerDegree  = p.Results.micronsPerDegree;
            obj.coneDarkNoiseRate = [0 0 0];
            obj.noiseFlag   = p.Results.noiseFlag;
            obj.emPositions = p.Results.emPositions;
            obj.useParfor   = p.Results.useParfor;
            
            % Construct outersgement if not passed.
            if (isempty(p.Results.os))
                eccentricityMeters = norm(p.Results.center);
                eccentricityDegs = eccentricityMeters / ...
                    (obj.micronsPerDegree / 1e6);
                obj.os = osLinear('eccentricity', eccentricityDegs);
            else
                obj.os = p.Results.os;
            end

            % Set the cone spacing and aperture given eccentricity & angle.
            %
            % Units returned are meters.
            %
            % passing p.Unmatched first ensures that we override it with
            % the expected units of meters and degrees used here,
            % independent of what the user passes.  Since the user can't
            % currently pass either eccentricity or angle, that is OK.
            % Because center is used in other places in coneMosaic, we
            % can't just override it with the passed 'eccentricity' and
            % 'angle' key value pairs.
            ecc = sqrt(sum(obj.center .^ 2));
            ang = atan2d(obj.center(2), obj.center(1));
            [spacing, aperture] = coneSizeReadData( p.Unmatched, ...
                'eccentricity', ecc, ...
                'eccentricityUnits', p.Results.eccentricityunits, ...
                'angle', ang, 'angleUnits', 'deg', ...
                'whichEye', obj.whichEye, ...
                'useParfor', obj.useParfor);

            obj.pigment.pdWidth  = aperture;
            obj.pigment.pdHeight = aperture;
            obj.pigment.height   = spacing;
            obj.pigment.width    = spacing;

            % See description of this parameter on the wiki page at
            % <https://github.com/isetbio/isetbio/wiki/
            %   Cone-mosaic#note-on-the-patternsamplesize-variable>
            obj.patternSampleSize = ...
                [obj.pigment.width, obj.pigment.height];

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
            
            % Return parameters if asked
            if nargout > 1
                % Return a struct with the parameters from this object.
                % If you want the default use
                %
                %   cmParams = coneMosaicRectP;
                
                cmParams.name    = obj.name;
                cmParams.pigment = obj.pigment; 
                cmParams.macular = obj.macular; 
                cmParams.os      = obj.os;
                cmParams.center  = obj.center;
                cmParams.wave    = obj.wave;
                cmParams.pattern = obj.pattern;
                cmParams.spatialDensity = obj.spatialDensity;
                cmParams.size    = obj.size;
                cmParams.integrationTime  = obj.integrationTime;
                cmParams.micronsPerDegree = obj.micronsPerDegree;
                cmParams.emPositions      = obj.emPositions;
                cmParams.apertureBlur     = obj.apertureBlur;
                cmParams.noiseFlag        = obj.noiseFlag; 
                cmParams.eccentricityunits = p.Results.eccentricityunits; 
                cmParams.whichEye = obj.whichEye; 
                cmParams.useParFor = obj.useParfor;            
            end
            
        end

        %% Get methods for dependent variables
        % <http://www.mathworks.com/help/matlab/matlab_oop/
        %   specifying-methods-and-functions.html#bu4wzba>
        % All functions that use dots in their names must be defined in the
        % classdef file, including:
        %   -  Converter methods that must use the package name as part of
        %      the class name because the class is contained in packages
        %   -  Property set and get access methods

        function val = get.wave(obj)
            % Retrieve the cone mosaic object's wave value
            %
            % Syntax:
            %   val = get.wave(obj)
            %
            % Description:
            %    Retrieve the wave value of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The wave value(s)
            %
            % Optional key/value pairs:
            %    None.
            %
            val = obj.pigment.wave;
        end

        function val = get.rows(obj)
            % Retrieve the cone mosaic object's rows value
            %
            % Syntax:
            %   val = get.rows(obj)
            %
            % Description:
            %    Retrieve the number of rows in the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The number of rows
            %
            % Optional key/value pairs:
            %    None.
            %
            val = size(obj.pattern, 1);
        end

        function val = get.cols(obj)
            % Retrieve the cone mosaic object's columns value
            %
            % Syntax:
            %   val = get.cols(obj)
            %
            % Description:
            %    Retrieve the number of columns in the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The number of columns
            %
            % Optional key/value pairs:
            %    None.
            %
            val = size(obj.pattern, 2);
        end

        function val = get.mosaicSize(obj)
            % Retrieve the cone mosaic object's mosaic size
            %
            % Syntax:
            %   val = get.mosaicSize(obj)
            %
            % Description:
            %    Retrieve the mosaic size of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The mosaic size in a 1 x 2 vector.
            %
            % Optional key/value pairs:
            %    None.
            %
            val = [obj.rows obj.cols];
        end

        function val = get.patternSupport(obj)
            % Retrieve the cone mosaic object's pattern support
            %
            % Syntax:
            %   val = get.patternSupport(obj)
            %
            % Description:
            %    Retrieve the pattern support of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The wave value(s)
            %
            % Optional key/value pairs:
            %    None.
            %
            x = (1:obj.cols) * obj.patternSampleSize(1);
            x = x - mean(x);
            y = (1:obj.rows) * obj.patternSampleSize(2);
            y = y - mean(y);
            [xx, yy] = meshgrid(x, y);
            val(:, :, 1) = xx;
            val(:, :, 2) = yy;
        end

        function val = size(obj, varargin)
            % Return size of cone mosaic patch in meters (height, width)
            %
            % Syntax:
            %   val = size(obj, [varargin])
            %
            % Description:
            %    Retrieve the size (in [height, width] format) of the cone
            %    mosaic patch in meters.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The wave value(s)
            %
            % Optional key/value pairs:
            %    None yet.
            %
            % Notes:
            %    * TODO: Should allow varargin to specify unit scale, e.g.
            %      cMosaic.size('unit', 'microns')
            %
            val = [obj.height, obj.width];
        end

        function val = get.width(obj)
            % Retrieve the cone mosaic object's width
            %
            % Syntax:
            %   val = get.width(obj)
            %
            % Description:
            %    Retrieve the width of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The width.
            %
            % Optional key/value pairs:
            %    None.
            %
            val = obj.cols * obj.patternSampleSize(1);
        end

        function val = get.height(obj)
            % Retrieve the cone mosaic object's height
            %
            % Syntax:
            %   val = get.height(obj)
            %
            % Description:
            %    Retrieve the height of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The height
            %
            % Optional key/value pairs:
            %    None.
            %
            val = obj.rows * obj.patternSampleSize(2);
        end

        function val = get.fov(obj)
            % Retrieve the cone mosaic object's field of view
            %
            % Syntax:
            %   val = get.fov(obj)
            %
            % Description:
            %    Retrieve the field of view of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The field of view
            %
            % Optional key/value pairs:
            %    None.
            %
            focalLengthMeters = obj.micronsPerDegree / (2 * 1e6) ...
                / tand(0.5);
            val = 2 * ...
                atand([obj.width obj.height] / 2 / focalLengthMeters);
        end

        function val = get.coverage(obj)
            % Retrieve the cone mosaic object's coverage factor
            %
            % Syntax:
            %   val = get.coverage(obj)
            %
            % Description:
            %    Compute the cone mosaic's cone coverage factor
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The mosaic's coverage factor.
            %
            % Optional key/value pairs:
            %    None.
            %
            if (isa(obj, 'coneMosaicHex'))
                [~, geometricCoverage] = obj.retinalCoverage();
                val = geometricCoverage;
            else
                mosaicAreaInMeters = prod(obj.size);
                conesNum = numel(find(obj.pattern>1));
                coneAreaInMeters = obj.pigment.area;
                val = conesNum * coneAreaInMeters / mosaicAreaInMeters;
            end
        end
        
        function val = get.innerSegmentCoverage(obj)
            % Retrieve the cone mosaic object's coverage factor
            %
            % Syntax:
            %   val = get.innerSegmentCoverage(obj)
            %
            % Description:
            %    Compute the cone mosaic's inner segment coverage factor
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The mosaic's inner segment coverage factor.
            %
            % Optional key/value pairs:
            %    None.
            %
            if (isa(obj, 'coneMosaicHex'))
                [lightCollectingCoverage, ~] = obj.retinalCoverage();
                val = lightCollectingCoverage;
            else
                mosaicAreaInMeters = prod(obj.size);
                conesNum = numel(find(obj.pattern>1));
                coneAreaInMeters = obj.pigment.pdArea;
                val = conesNum * coneAreaInMeters / mosaicAreaInMeters;
            end
        end
        
        function val = get.coneLocs(obj)
            % Retrieve the cone mosaic object's cone locations (in meters)
            %
            % Syntax:
            %   val = get.coneLocs(obj)
            %
            % Description:
            %    Retrieve the cone locations of the cone mosaic object and
            %    return them in meters.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The cone locations in meters.
            %
            % Optional key/value pairs:
            %    None.
            %
            x = (1:obj.cols) * obj.pigment.width;
            x = x - mean(x) + obj.center(1);
            y = (1:obj.rows) * obj.pigment.height;
            y = y - mean(y) + obj.center(2);

            [X, Y] = meshgrid(x, y);
            val = [X(:) Y(:)];
        end

        function val = get.qe(obj)
            % Return effective absorptance including macular pigment
            %
            % Syntax:
            %   val = get.qe(obj)
            %
            % Description:
            %    Compute and return the effective absorptance of a cone
            %    mosaic object, including macular pigment.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - effective absorptance
            %
            % Optional key/value pairs:
            %    None.
            %
            
            % val = bsxfun(@times, obj.pigment.absorptance, ...
            %    obj.macular.transmittance) * ...
            %    diag(obj.pigment.peakEfficiency);
            val = bsxfun(@times, obj.pigment.quantalEfficiency, ...
                obj.macular.transmittance) ;
        end

        function val = get.spatialDensity(obj)
            % Retrieve the cone mosaic object's spatial density
            %
            % Syntax:
            %   val = get.spatialDensity(obj)
            %
            % Description:
            %    Retrieve the spatial density of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The spatial density
            %
            % Optional key/value pairs:
            %    None.
            %
            val = obj.spatialDensity_;
        end

        % Shouldn't have to do this any more, right?
        function val = get.current(obj)
            % Retrieve the cone mosaic object's current
            %
            % Syntax:
            %   val = get.current(obj)
            %
            % Description:
            %    Retrieve the current of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The current
            %
            % Optional key/value pairs:
            %    None.
            %
            if isempty(obj.current)
                val = [];
            else
                val = double(obj.current);
            end
        end

        function val = get.absorptions(obj)
            % Retrieve the cone mosaic object's absorptions
            %
            % Syntax:
            %   val = get.absorptions(obj)
            %
            % Description:
            %    Retrieve the absorptions of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The absorptions
            %
            % Optional key/value pairs:
            %    None.
            %
            val = double(obj.absorptions);
        end

        function val = get.tSamples(obj)
            % Return computed tSamples from 1st dim. of emPositions
            %
            % Syntax:
            %   val = get.tSamples(obj)
            %
            % Description:
            %    Retrieve the computed tSamples from the first dimension of
            %    the eye movement positions.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The t samples
            %
            % Optional key/value pairs:
            %    None.
            %
            if isempty(obj.emPositions), val = 1; return; end
            
            % The number of eye positions defines the number of samples
            val = size(obj.emPositions, 1);
        end

        function val = timeStep(obj)
            % Retrieve the cone mosaic object's current time step.
            %
            % Syntax:
            %   val = get.wave(obj)
            %
            % Description:
            %    Retrieve the current time step of the cone mosaic object
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The current time step
            %
            % Optional key/value pairs:
            %    None.
            %
            val = obj.os.timeStep;
        end

        %% Set methods for class properties
        function set.spatialDensity(obj, val)
            % Assign a value to the cone mosaic object's spatial density
            %
            % Syntax:
            %   set.spatialDensity(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's spatial density.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired spatial density
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            if all(obj.spatialDensity_(:) == val(:)), return; end
            obj.spatialDensity_ = val;
            [~, obj.pattern] = humanConeMosaic(obj.mosaicSize, ...
                val, obj.patternSampleSize(1));
            obj.clearData();
        end

        function set.wave(obj, val)
            % Assign a value to the cone mosaic object's wave parameter
            %
            % Syntax:
            %   set.wave(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's wave parameter. This will apply to both the
            %    pigment and macular aspects.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired wavelength(s)
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.pigment.wave = val(:);
            obj.macular.wave = val(:);
        end

        function set.absorptions(obj, val)
            % Assign a value to the cone mosaic object's absorptions
            %
            % Syntax:
            %   set.absorptions(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's absorptions.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired absorptions
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.absorptions = single(val);
        end

        function set.integrationTime(obj, val)
            % Assign a value to the cone mosaic object's integration time.
            %
            % Syntax:
            %   set.integrationTime(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's integrationTime.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired integration time
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.integrationTime = val;
        end

        function set.current(obj, val)
            % Assign a value to the cone mosaic object's current
            %
            % Syntax:
            %   set.current(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's current.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired current.
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.current = single(val);
        end

        function set.fov(obj, val)
            % Assign a value to the cone mosaic object's field of view
            %
            % Syntax:
            %   set.fov(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's field of view.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired field of view
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            if (isa(obj, 'coneMosaicHex'))
                error('Setting the ''fov'' property directly is not allowed for a coneMosaicHex object. The fov is settable only during instantiation.');
            else
                obj.setSizeToFOV(val);
            end
        end

        function set.mosaicSize(obj, val)
            % Assign a value to the cone mosaic object's mosaic size
            %
            % Syntax:
            %   set.mosaicSize(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's mosaic size.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired mosaic size
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.spatialDensity_, obj.patternSampleSize(1));
                obj.os.patchSize = obj.width;
                obj.clearData();
            end
        end

        function set.rows(obj, val)
            % Assign a value to the cone mosaic object's rows parameter
            %
            % Syntax:
            %   set.rows(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's rows parameter
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired rows value
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.mosaicSize = [val obj.cols];
        end

        function set.cols(obj, val)
            % Assign a value to the cone mosaic object's columns parameter
            %
            % Syntax:
            %   set.cols(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's columns parameter.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired columns parameter.
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            obj.mosaicSize = [obj.rows val];
        end

        function set.noiseFlag(obj, val)
            % Assign a value to the cone mosaic object's noise flag.
            %
            % Syntax:
            %   set.noiseFlag(obj, val)
            %
            % Description:
            %    Assign the provided value (val) to the cone mosaic
            %    object's noise flag after confirming it is on the list of
            %    possible valid flags.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired noise flag.
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            if ischar(val) && ...
                    (ismember(lower(val), coneMosaic.validNoiseFlags))
                obj.noiseFlag = val;
            else
                s = sprintf('%s ', coneMosaic.validNoiseFlags{:});
                error(['''%s'' is an invalid value for ' ...
                    'coneMosaic.noiseFlag. Choose one from: %s '], val, s);
            end
        end
    end

    methods (Access=public)
        % Declare the compute method
        [absorptions, current, varargout] = compute(obj, oi, varargin);

        % Declare the compute method for a sequence of optical images
        % viewed sequentially.
        [absorptions, absorptionsTimeAxis, varargout] = ...
            computeForOISequence(obj, oiSequence, varargin)

        function [demosaicedAbsorptionsMap, sRGB] = ...
                demosaicedIsomerizationMaps(obj, varargin)
            % Return demosaiced isomerization maps and corresponding sRGB
            %
            % Syntax:
            %   [demosaicedAbsorptionsMap, sRGB] = ...
            %       demosaicedIsomerizationMaps(obj, [varargin])
            %
            % Description:
            %    This method will return the demosaiced isomerization maps
            %    and the corresponding sRGB rendition.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %    val - The desired spatial density
            %
            % Outputs:
            %    None
            %
            % Optional key/value pairs:
            %    None.
            %
            [demosaicedAbsorptionsMap, sRGB] = obj.demosaicedResponses();
        end

        function demosaicedCurrentsMap = ...
                demosaicedPhotoCurrentMaps(obj, current)
            % Return the demosaiced photocurrent maps
            %
            % Syntax:
            %   demosaicedCurrentMap = demosaicedPhotoCurrentMaps(obj, ...
            %       current)
            %
            % Description:
            %    Return the demosaiced photocurrent maps
            %
            % Inputs:
            %    obj                   - The cone mosaic object
            %    current               - The current value
            %
            % Outputs:
            %    demosaicedCurrentsMap - The demosaiced photocurrent maps
            %
            % Optional key/value pairs:
            %    None.
            %
            demosaicedCurrentsMap = obj.demosaicedResponses(current);
        end

        % Demosaicing method
        varargout = demosaicedResponses(obj, varargin);

        % Method that low-passes a mosaic's response using different space
        % constants for each of the L, M, and S-cone submosaic.
        [lowPassedResponse, Lmap, Mmap, Smap] = ...
            lowPassActivationMap(response, theMosaic, spaceConstants);

        % Declare the computeCurrent method
        [interpFilters, meanCur] = computeCurrent(obj, varargin);

        % Declare method for computing retinal coverage
        [apertureCoverage, geometricCoverage] = retinalCoverage(obj);
        
        % Determine whether we should correct absorptions with eccentricity
        correctAbsorptionsForEccentricity = ...
            shouldCorrectAbsorptionsWithEccentricity(obj);
        
        function val = timeAxis(obj)
            % Return the cone mosaic time axis
            %
            % Syntax:
            %   val = timeAxis(obj)
            %
            % Description:
            %    Return the cone mosaic's time axis
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - the time axis value
            %
            % Optional key/value pairs:
            %    None.
            %
            val = (0:1:(obj.tSamples - 1)) * obj.integrationTime;
        end

        function val = interpFilterTimeAxis(obj)
            % Return the interpolation filter's time axis
            %
            % Syntax:
            %   val = interpFilterTimeAxis(obj)
            %
            % Description:
            %    Return the time axis for the interpolation filter of the
            %    cone mosaic object.
            %
            % Inputs:
            %    obj - The cone mosaic object
            %
            % Outputs:
            %    val - The interpolation filter's time axis
            %
            % Optional key/value pairs:
            %    None.
            %
            osTimeAxis = obj.os.timeAxis;
            if (isempty(osTimeAxis))
                val = [];
            else
                dT = obj.integrationTime;
                val = (osTimeAxis(1):dT:osTimeAxis(end));
            end
        end

    end

    methods (Static)
        [noisyImage, theNoise] = photonNoise(absorptions, varargin);
        resampledAbsorptionsSequence = tResample(absorptionsSequence, ...
            pattern, originalTimeAxis, resampledTimeAxis);
        validateClippingRect(clippingRect);
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