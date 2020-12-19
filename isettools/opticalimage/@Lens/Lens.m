classdef Lens < handle
% Class for human lens pigment properties
%
% Syntax:
%   obj = Lens;
%
% Description:
%    Creates a lens object to describe absorbance properties of the
%    human lens.
%
%    Our understanding of the terminology and best conventions for
%    describing this type of data has evolved over time, and the property
%    names differ somewhat from the conventions we would adopt today.
%    Changing tne names in the code will produce backwards compatibility
%    issues, however, so we have done our best to comment and explain here.
%
%    obj.unitDensity is the lens absorbance spectrum. It is called unit
%    density for historical reasons.  We multiply this by a scalar
%    (obj.density), described just below before using this to compute. It
%    might be better to call the product the absorbance spectrum, how this
%    nomenclature works seems to vary in the literature. Note that in other
%    similar classes (e.g., Macular), the absorbance in obj.unitDensity is
%    normalized to a peak of 1, but not here.
%
%    obj.density is a scalar applied to the absorbance stored in the
%    unitDensity property. Because this number multiplies the unit density
%    property when either is used, the two quantities need to be
%    interpretted together to give the actual absorbance. Often one follows
%    the convention that unit density is normalized to a max of 1, but that
%    is not done here, and indeed given the data we read in for unit
%    density and the fact that we don't scale it, the default value for
%    obj.density in this routine is 1. This is in contrast to the way the
%    very similar Macular and photoPigment objects are coded, where the
%    unit density is scaled to a maximum of 1 and the obj.density parameter
%    is interpretted as the peak optical density.
%
%    obj.absorptance is the absorptance spectrum.
%
%    obj.transmittance is the transmission spectrum.
%
%    The absorbance data that drive this routine are stored on wavelength
%    support in property wave_ in property unitDensity_.  Typically this is
%    set to a large wavelength support and then interpolated onto the
%    support in propety wave.  You can set unitDensity_ after the object is
%    instantiated, but you can't change wave_.  When you set unitDensity,
%    it should be on wavelength support wave, and it is splined onto
%    the wavelength support in wave_ before being stored in unitDensity_.
%    Note that this design does not prevent you from setting unitDensity on
%    wavelength support very different from that being used to store the
%    data, which could lead to extrapolation errors.  To avoid this, if you
%    want to use custom data, you may be better off creating the object
%    with the desired data on the wavelength support you intend to use.
%    That said, the default values are read in and stored on wave_ support
%    of 390:830 at 1 nm spacing, which is good for most applications.
%
%    Useful formulae:
%       Absorbance spectra are normalized to a peak value of 1. In this
%       routine, for historical reasons, absorbance is called unitDensity.
%
%       Absorptance spectra are the proportion of quanta actually absorbed.
%       This is the term used in this routine.
%       
%       Equation: absorptance = 1 - 10.^(-opticalDensity * absorbance).  In
%       this routine, again for historical reasons, opticalDensity is just
%       called density.  In the literature, this is sometimes called peak
%       optical density.
%
%       Transmittance is 1-absorptance, the amount of light that passes
%       through the pigment.  Alternately, you can compute
%       transmittance = 10.^(-opticalDensity * absorbance) and
%       absorptance as = 1-transmittance.
%
%    The default macular density absorbance was obtained from an old Stockman
%    site, but should match that at the new Stockman site (cvrl.org) and
%    those in the Psychtoolbox.
%
% Inputs:
%    None required.
%
% Outputs:
%    lens        - Object. The lens object.
%
% Optional key/value pairs:
%    name        - String. The name for this object, default 'human lens'.
%    wave        - Vector. The wavelength samples.
%    density     - Numeric. The lens pigment density, default 1.
%    unitDensity - Numeric. The unit density. Default [], which causes the
%                  routine to reads the densities from lensDensity.mat and
%                  sets wave_ to the native wavelenght support of that
%                  file, 390:830. If passed, should be on same wavelength
%                  support as wave, and the wave_ property is set to wave.
%
% See Also:
%   opticsGet, opticsSet, Macular, photoPigment
%

% History:
%    xx/xx/13  HJ/BW  ISETBIO Team 2013.
%    03/05/18  jnm    Formatting
%    07/02/19  JNM    Formatting update
%    12/13/20  dhb    Comments and code cleaning.
%    12/18/20  dhb    More comments.

% Examples:
%{
    lens = Lens;
%}
%{
    % Specify an observer with lens absorbance 20% higher than the default
    % value.
    thisLens = Lens('wave', 400:10:700, 'density', 1.2, 'name', 'my aging lens');
%}
%{
    % Show that abosrbtance and transmission sum to 1.
    lens = Lens;
    t = lens.transmittance;
    a = lens.absorptance;
    vcNewGraphWin;
    plot(lens.wave, t, 'r-', lens.wave, a, 'g-', lens.wave, a + t, 'k--');
    xlabel('Wave (nm)');
    ylabel('Fraction');
    legend({'transmittance', 'absorptance', 'sum'});
%}

properties  % public properties
    % name - Name of this particular lens object
    name;

    % density - pigment density
    density;

    % wave - wavelength samples in nm
    wave;
end

properties (SetAccess=private)
    % wave_ - internal wavelength samples
    wave_;

    % unitDensity_ - unit density absorbance sampled with wave_
    unitDensity_;
end

properties (Dependent)
    % unitDensity - spectral absorbance with unit pigment density
    unitDensity;

    % spectralDensity - unitDensity scaled by obj.density
    spectralDensity;

    % transmittance - proportion of quanta transmitted
    transmittance;

    % absorptance - proportion of quanta absorbed
    absorptance;
end

methods  % public methods
    % constructor
    function obj = Lens(varargin)
    % Lens constructor. See above for more information.
    
    % parse input
    p = inputParser;
    p.addParameter('wave', 400:10:700, @isnumeric);
    p.addParameter('density', 1, @isscalar);
    p.addParameter('unitDensity', [], @isnumeric);
    p.addParameter('name', 'human lens', @isstr);
    p.parse(varargin{:});

    % set properties
    obj.name = p.Results.name;
    obj.wave = p.Results.wave(:);
    obj.density = p.Results.density;

    if isempty(p.Results.unitDensity)
        obj.wave_ = (390:830)';
        obj.unitDensity_ = ieReadSpectra('lensDensity.mat', ...
            obj.wave_);
    else
        obj.wave_ = p.Results.wave;
        obj.unitDensity_ = p.Results.unitDensity;
    end
    end

    % get methods for dependent variables
    function val = get.unitDensity(obj)
    % Retrieve Lens unit density
    %
    % Syntax:
    %   val = get.unitDensity(obj)
    %
    % Description:
    %    Retrieve the Lens object's unit density
    %
    % Inputs:
    %    obj - Object. The lens object.
    %
    % Outputs:
    %    val - Numeric. The unit density.
    %
    % Optional key/value pairs:
    %    None.
    %

    % interpolate for wavelength samples
    val = interp1(obj.wave_, obj.unitDensity_, obj.wave, 'pchip');
    val = max(val, 0);
    end

    function val = get.spectralDensity(obj)
    % Retrieve the Lens object spectral density
    %
    % Syntax:
    %   val = get.spectralDensity(obj)
    %
    % Description:
    %    Retrieve the Lens object's spectral density
    %
    % Inputs:
    %    obj - Object. The lens object.
    %
    % Outputs:
    %    val - Numeric. The spectral density.
    %
    % Optional key/value pairs:
    %    None.
    %

    % compute scaled absorbance
    val = obj.unitDensity * obj.density;
    end

    function val = get.transmittance(obj)
    % Compute the proportion of quanta transmitted
    %
    % Syntax:
    %   val = get.transmittance(obj)
    %
    % Description:
    %    Retrieve the Lens object's transmittance
    %
    % Inputs:
    %    obj - Object. The lens object
    %
    % Outputs:
    %    val - The transmittance
    %
    % Optional key/value pairs:
    %    None.
    %

    % Can interpolate
    val = 10 .^ (-obj.spectralDensity);
    end

    function val = get.absorptance(obj)
    % compute proportion of quanta absorbed
    %
    % Syntax:
    %   val = get.absorptance(obj)
    %
    % Description:
    %    Retrieve the Lens object's proportion of quanta absorbed
    %
    % Inputs:
    %    obj - The lens object
    %
    % Outputs:
    %    val - The absorptance
    %
    % Optional key/value pairs:
    %    None.
    %

    val = 1 - obj.transmittance;
    end

    function val = get.density(obj)
    % Retrieve the Lens object's density
    %
    % Syntax:
    %   val = get.density(obj)
    %
    % Description:
    %    Retrieve the Lens object's density
    %
    % Inputs:
    %    obj - The lens object
    %
    % Outputs:
    %    val - The density
    %
    % Optional key/value pairs:
    %    None.
    %

    val = obj.density;
    end

    function set(obj, param, val)
    % Assign the provided value to the specified parameter
    %
    % Syntax:
    %   set(obj, param, val)
    %
    % Description:
    %    Using the provided values, assign or change the value of the
    %    specified parameter to that of the provided value.
    %
    % Inputs:
    %    obj   - The lens object
    %    param - String. The parameter you wish to modify. Options
    %            include the following:
    %      'name'        - String. The Lens object name.
    %      'wave'        - Numerical. The wavelength vector.
    %      'density'     - Scalar value. The lens pigment density.
    %
    %     val   - The value to assign
    %
    % Outputs:
    %    None.
    %
    % Optional key/value pairs:
    %    None.
    %

    p = inputParser;
    p.KeepUnmatched = true;
    p.addRequired('param', @isstr);
    p.addRequired('val');

    % set parameter value
    switch ieParamFormat(param)
        case 'name'
            assert(ischar(val), 'Name should be a string');
            obj.name = val;
        case {'wave', 'wavelength'}
            assert(isvector(val), 'wave should be vector');
            obj.wave = val(:);
        case {'absorbance', 'unitdensity'}
            assert(length(val) == length(obj.wave), ...
                'Val should have same length as lens wavelength');
            obj.unitDensity_ = interp1(...
                obj.wave, val, obj.wave_, 'pchip');
            obj.unitDensity_ = max(obj.unitDensity_, 0);
        case 'density'
            assert(isscalar(val), 'val should be scalar');
            obj.density = val;
        otherwise
            error('Unknown parameter %s\n', param);
    end
    end

end
end