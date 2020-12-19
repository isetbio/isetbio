classdef photoPigment < hiddenHandle
% Class for single cone photopigment and related properties
%
% Syntax:
%	pigment = photoPigment;
%
% Description:
%    This class contains properties for the photopigment absorption
%    properties of a single cone cell. For the full cone mosaic, see the
%    coneMosaic and coneMosaicHex classes.
%
%    Our understanding of the terminology and best conventions for
%    describing this type of data has evolved over time, and the property
%    names differ somewhat from the conventions we would adopt today.
%    Changing tne names in the code will produce backwards compatibility
%    issues, however, so we have done our best to comment and explain here.
%
%    Most of the terms represented here are descriptions of the
%    photopigment itself. In addition, there are a few terms the
%    capture the effective optical size of the photopigment absorption.
%
%    obj.absorbance is the absorbance spectra of the L, M, and S cones, each
%    normalized to a peak value of 1. Values are in the columns, with a separate
%    column for L, M and S. The column-wise arrangement applies to this and to 
%    the other L, M, and S spectra below. This quantity is called obj.unitDensity
%    in the Lens and Macular objects, and in the Lens object it is not normalized.
%    The normalization to peak of 1 is just a convention, the quantity that
%    matters is the product obj.opticalDensity*ojb.absorbance, and there are
%    those who might call that product the absorbance.
%
%    obj.opticalDensity is the peak optical density (sometimes just called optical
%    density). The interpretation as peak optical density depends on the
%    convention followed here of normlizing obj.absorbance to a peak of 1.
%    There are three entries to this vector, one each for the L, M, and S
%    cones.
%
%    obj.peakEfficiency is the probability that a photopigment absorption leads
%    to an isomerization.  There are three entries to this vector, one each for the L,
%    M, and S cones.  This property is unfortunatley named, as is it isn't the
%    peak of anything.
%
%    obj.absorptance is the absorptance spectrum.  This tells us the
%    probability that a photon of a given wavelength is absorbed as it
%    passes through a layer of photopigment with total absorbance given
%    by obj.opticalDensity*ojb.absorbance.
%
%    obj.quantalEfficiency. These are the actual quantal efficiences with
%    which an incident photon causes an isomerization.  Obtained by
%    multiplying obj.absorptance by the quantal efficiency
%    obj.peakEfficiency.  These are at the cone, and do not take into
%    account effect of lens or macular pigment, nor of cone collecting
%    area.  This is true of the fundamentals below.
% 
%    obj.quantaFundamentals are nomalized (each to a peak of 1) fundamentals
%    of the L, M, and S cones, in quantal units. These really shouldn't be
%    used for anything other than perhaps making a plot of the relative
%    shapes of the photopigment action spectra expressed in quantal units.
%    The term cone fundamentals usually is takent to mean cone
%    sensitivities expressed relative to light entering the eye ("at the
%    cornea"), which these are not.  They are also not in any useful units,
%    because of the normalization.
%
%    obj.energyFundamentals. Same as obj.quantaFundamentals, but in energy
%    units. As with obj.quantaFundamentals, these are not in useful units
%    because they begin with the normalized obj.quantaFundamentals.  So,
%    possibly of interest for seeing the shape of the function, but not
%    anything one should be encouraged to use.
%
%    Useful formulae:
%       Absorbance spectra here are normalized to a peak value of 1, and
%       then scaled by optical density to get the not normalized
%       absorbance.
%
%       Absorptance spectra are the proportion of quanta actually absorbed.
%       This is the term used in this routine.
%       
%       Equation: absorptance = 1 - 10.^(-opticalDensity * absorbance).  In
%       this routine, again for historical reasons, opticalDensity is just
%       called density.  In the literature, this is sometimes called peak
%       optical density.
%
%       The absorptance
%
%    The absorbance data that drive this routine are stored on wavelength
%    support in property wave_ in property absorbance.  Typically wave_ is
%    set to a large wavelength support and then interpolated onto the
%    support in propety wave.  You can set absorbance after the object is
%    instantiated, but you can't change wave_.  When you set absorbance, it
%    should be on wavelength support wave, and it is splined onto the
%    wavelength support in wave_ before being stored in obj.absorbance.
%    Note that this design does not prevent you from setting absorbance on
%    wavelength support very different from that being previously used to
%    store the data, which could lead to extrapolation errors.  To avoid
%    this, if you want to use custom data, you may be better off creating
%    the object with the desired data on the wavelength support you intend
%    to use. That said, the default values are read in and stored on wave_
%    support of 390:830 at 1 nm spacing, which is good for most
%    applications.
%
% Input:
%	 None required.
%
% Output:
%    pigment          - The created photoPigment object.
%   
% Optional key/value pairs:
%	 'wave'           - Vector of wavelengths in nm (400:10:31).
%    'opticalDensity' - Three vector of peak optical densities for L, M and
%                       S cone photopigment. Default: [0.5 0.5 0.4].
%    'absorbance'     - L, M and S cone absorbance spectra. Default
%                       empty, in which case these are obtained through
%                       routine coneAbsorbanceReadData.
%    'peakEfficiency' - Quantal efficiency for isomerizations for
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
%    t_conePhotoPigment, cPhotoPigment, coneMosaic, Macular, Lens
%

% History:
%    xx/xx/16  HJ   ISETBIO Team, 2016
%    02/15/18  jnm  Formatting
%    12/18/20  dhb  Comments.  Add quantalEfficiency property.

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

    % quantaFundamentals - normalized cone absorptance. Because of the
    % normalization, not useful for actual calculation, only for examining
    % the relative shape.
    quantaFundamentals;

    % energyFundamentals - normalized cone absorptance converrted
    % to energy units from the normalized obj.quantaFundamentals.
    % Useful only for shape, because of the normalization.
    energyFundamentals;
    
    % quantalEfficiency - actual probability of isomerization in real
    % quantal units. Gives the probability that an incident photon
    % isomerizes photopigment. These are actually useful.  Does not take
    % into account effect of inert pigments (lens, macular pigment).
    quantalEfficiency;

    % area - The area of the object. Calculated by width * height
    area;

    % pdArea - The pdArea of the object. Calculated by pdWidth * pdHeight
    pdArea;

    % gapWidth - the width of the gap. Calculated by width - pdWidth
    gapWidth;

    % gapHeight - The height of the gap. Calculated by height - pdHeight
    gapHeight;
end

properties (SetAccess=private)
    % wave_ - The internal wavelength samples
    wave_;

    % absorbance_ - The absorbance data sampled at wave_
    absorbance_;
end

methods  % public methods
    % constructor
    function obj = photoPigment(varargin)
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

    function val = get.quantalEfficiency(obj) % compute absorptance
        % Retrieve photo pigment object's isomerization efficiency.
        %
        % Syntax:
        %   obj = get.quantalEfficiency(obj)
        %
        % Description:
        %    Retrieve the quantal efficiency from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The LMS quantal efficiencies for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.absorptance*diag(obj.peakEfficiency);
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
    function set.absorbance(obj, val)
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
        %          sampling (wave, not wave_). The passed values are
        %          clipped to range 0 to 1.
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    None.
        %
        
        % Set
        obj.absorbance_ = interp1(obj.wave, val, obj.wave_, ...
            'linear', 'extrap');
        obj.absorbance_ = ieClip(obj.absorbance_, 0, 1);
    end
end

methods (Static)
    % When we have them, they go here
end
end