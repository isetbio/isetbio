classdef photoPigment < hiddenHandle
% Class for single cone photopigment and related properties
%
% Syntax:
%	pigment = photoPigment;
%
% Description:
%    This class contains properties for the photopigment absorption
%    properties of a single cone cell. 
%
%    For the full cone mosaic, see the coneMosaic class
%
%    Most of the terms represented here are descriptions of the
%    photopigment itself. In addition, there are a few terms the
%    capture the effective optical size of the photopigment absorption.
%
%    Default parameters are determined by underlying routines that get
%    the required data types. Unmatched key/value pairs passed to
%    photoPigment are passed on to the underlying routines and can be
%    used to adjust the parameters obtained. See help for each routine
%    for what the available key/value pairs are.
%
%         absorbance    coneAbsorbanceReadData
%
% Input:
%	 None required.
%
% Output:
%    pigment          - The created photoPigment object.
%   
% Optional key/value pairs:
%	 'wave'           - Vector of wavelengths in nm (400:10:31).
%    'opticalDensity' - Three vector of optical densities for L, M and
%                       S cone photopigment. Default: [0.5 0.5 0.4].
%    'absorbance'     - L, M and S cone absorbance spectra. Default
%                       empty, in which case these are obtained through
%                       routine coneAbsorbanceReadData.
%    'peakEfficiency' - Peak quantal efficiency for isomerizations for
%                       L, M and S cones. Default [2 2 2]/3.
%    'width'          - Cone width (including gap between cones) in
%                       meters. Default 2e-6.
%    'height'         - Cone height (including gap between cones) in
%                       meters. Default 2e-6.
%    'pdWidth'        - Collecting area width in meters (default 2e-6)
%    'pdHeight'       - Collecting area height in meters (default 2e-6)
%
% Notes:
%    * [NOTE: DHB - Need to explain about width and height, pdWidth and
%      pdHeight and how these are used. Perhaps even simplify code not
%      to have both.]
%
% See Also:
%    t_conePhotoPigment, coneMosaic, Macular, lens
%

% History:
%    xx/xx/16  HJ   ISETBIO Team, 2016
%    02/15/18  jnm  Formatting

properties  % public properties
    % opticalDensity - photopigment optical densities for L, M, S
    opticalDensity;

    % peakEfficiency - peak absorptance efficiency
    peakEfficiency;

    % width - cone width (including gap) in meters
    width;

    % height - cone height (including gap) in meters
    height;

    % pdWidth - photodetector width in meters
    pdWidth;

    % pdHeight - photodetector height in meters
    pdHeight;
end

properties (SetObservable, AbortSet)
    % wave - wavelength samples
    wave;
end

properties (Dependent)
    % absorbance - spectral absorbance of the cones
    absorbance;

    % absorptance - cone absorptance without ocular media
    absorptance;

    % quantaFundamentals - normalized cone absorptance
    quantaFundamentals;

    % energyFundamentals - normalized cone absorption in energy units
    energyFundamentals;

    % area - The area of the object. Calculated by width * height
    area;

    % pdArea - The pdArea of the object. Calculated by pdWidth * pdHeight
    pdArea;

    % gapWidth - the width of the gap. Calculated by width - pdWidth
    gapWidth;

    % gapHeight - The height of the gap. Calculated by height - pdHeight
    gapHeight;
end

properties(Access = private)  % private properties
    % wave_ - The internal wavelength samples
    wave_;

    % absorbance_ - The absorbance data sampled at wave_
    absorbance_;
end

methods  % public methods
    % constructor
    function obj = photoPigment(varargin)
        % Initialize defaults for photoPigments parameters
        %
        % Syntax:
        %   obj = photoPigment([varargin]);
        %
        % Description:
        %    Initialize the default values for the public properties: wave
        %    (400:10:700), opticalDensity ([.5 .5 .4]), absorbance ([]),
        %    peakEfficiency ([2 2 2]/3), width (2e-6), height (2e-6),
        %    pdWidth (2e-6), and pdHeight (2e-6). And then for the
        %    dependent and private object properties.
        %
        % Inputs:
        %    None required.
        %
        % Outputs:
        %    obj - The created photo pigment object
        %
        % Optional key/value pairs:
        %    None.
        %
        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
        p.addParameter('absorbance', [], @isnumeric);
        p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
        p.addParameter('width', 2e-6, @isnumeric);
        p.addParameter('height', 2e-6, @isnumeric);
        p.addParameter('pdWidth', 2e-6, @isnumeric);
        p.addParameter('pdHeight', 2e-6, @isnumeric);

        p.parse(varargin{:});

        % set object properties
        obj.wave = p.Results.wave(:);
        obj.opticalDensity = p.Results.opticalDensity(:);
        obj.peakEfficiency = p.Results.peakEfficiency(:);

        obj.width = p.Results.width;
        obj.height = p.Results.height;
        obj.pdWidth = p.Results.pdWidth;
        obj.pdHeight = p.Results.pdHeight;

        % If absorbance is not specified, we obtain it using the defaults
        % of coneAbsorbanceReadData. 
        if isempty(p.Results.absorbance)
            obj.wave_ = (390:830)';
            obj.absorbance_ = coneAbsorbanceReadData(p.Unmatched, ...
                'wave', obj.wave_);
        else
            obj.wave_ = p.Results.wave;
            obj.absorbance_ = p.Results.absorbance;
        end
    end

    % get method for dependent variable
    function val = get.absorbance(obj) % interpolate for absorbance
        % Retrieve photo pigment object's absorbance value
        %
        % Syntax:
        %   obj = get.absorbance(obj)
        %
        % Description:
        %    Retrieve the absorbance from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The absorbance value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = interp1(obj.wave_, obj.absorbance_, obj.wave, ...
            'linear', 'extrap');
        val = ieClip(val, 0, 1);
    end

    function val = get.absorptance(obj) % compute absorptance
        % Retrieve photo pigment object's absorptance value
        %
        % Syntax:
        %   obj = get.absorptance(obj)
        %
        % Description:
        %    Retrieve the absorptance from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The absorptance value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = 1 - 10 .^ (-obj.absorbance * diag(obj.opticalDensity));
    end

    function val = get.quantaFundamentals(obj)
        % compute and return quanta fundamentals
        %
        % Syntax:
        %   obj = get.quantaFundamentals(obj)
        %
        % Description:
        %    Compute and return the quanta fundamentals for the photo
        %    pigment object obj.
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The quanta fundamentals for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = bsxfun(@rdivide, obj.absorptance, max(obj.absorptance));
    end

    function val = get.energyFundamentals(obj)
        % Retrieve photo pigment object's energy fundamentals
        %
        % Syntax:
        %   obj = get.energyFundamentals(obj)
        %
        % Description:
        %    Retrieve the energy fundamentals from the photoPigment object
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The energy fundamentls for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        h = vcConstants('planck');
        c = vcConstants('speed of light');
        val = 1e-9 * bsxfun(@times, obj.quantaFundamentals / h / c, ...
            obj.wave);
    end

    function val = get.area(obj)
        % Retrieve photo pigment object's area
        %
        % Syntax:
        %   obj = get.area(obj)
        %
        % Description:
        %    Retrieve the area from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.width * obj.height;
    end

    function val = get.gapWidth(obj)
        % Retrieve photo pigment object's gap width value
        %
        % Syntax:
        %   obj = get.gapWidth(obj)
        %
        % Description:
        %    Retrieve the gap width from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The gap width value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.width - obj.pdWidth;
    end

    function val = get.gapHeight(obj)
        % Retrieve photo pigment object's gap height value
        %
        % Syntax:
        %   obj = get.gapHeight(obj)
        %
        % Description:
        %    Retrieve the gap height from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The gap height value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.height - obj.pdHeight;
    end

    function val = get.pdArea(obj)
        % Retrieve photo pigment object's pd area value
        %
        % Syntax:
        %   obj = get.pdArea(obj)
        %
        % Description:
        %    Retrieve the pd area from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The pd area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.pdWidth * obj.pdHeight;
    end

    % set method for dependent variable
    function set.absorbance(obj, val, varargin)
        % Set the photo pigment object's absorbance value
        %
        % Syntax:
        %   obj = set.absorbance(obj, val)
        %
        % Description:
        %    Set the photo pigment's absorbance spectrum.  Setting this
        %    overrides 
        %
        % Inputs:
        %    obj - The photoPigment object
        %    val - The absorbance value to set.  Should be on
        %          same wavelength spacing as the object's wavelength
        %          sampling (wave, not wave_).  Could add a key value pair
        %          to allow passing of the wavelength support being passed.
        %          The passed values are clipped to range 0 to 1.
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    'wave'     - Wavelength support of passed absorbance.  If not
        %                 passed, this is assumed to be the objects 'wave' 
        %                 support, and the absorbance is interpolated to the object's
        %                 'wave_' support.  This is set when the object is
        %                 created, and by default is 390:830. If it is
        %                 passed, then the object's 'wave_' support is
        %                 changed to the passed value.
        %
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', [], @isnumeric);
        p.parse(varargin{:});
        
        if (isempty(p.Results.wave))
            obj.absorbance_ = interp1(obj.wave, val, obj.wave_, ...
                'linear', 'extrap');
        else
            obj.absorbance = val;
            obj.wave_ = p.Results.wave;
        end
        
        % Clip into reasonable range
        obj.absorbance_ = ieClip(obj.absorbance_, 0, 1);
    end
end

methods (Static)
    % When we have them, they go here
end
end