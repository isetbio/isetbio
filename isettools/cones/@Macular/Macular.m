classdef Macular < hiddenHandle
% Class for macular pigment properties
%
% Syntax:
%   obj = Macular;
%
% Description:
%    The human retina contains a pigment that covers the central (macular)
%    region. This macular pigment passes certain wavelengths of light more
%    than others. The pigment varies in density from central vision, where
%    it is highest, to increasingly peripheral vision.
%
%    This class manages several measures of the macular pigment wavelength
%    properties as a function of macular pigment density. Our understanding
%    of the terminology and best conventions for describing this type of
%    data has evolved over time, and the property names differ somewhat
%    from the conventions we would adopt today. Changing tne names in the
%    code will produce backwards compatibility issues, however, so we have
%    done our best to comment and explain here.
%
%    obj.unitDensity is the absorbance spectrum, normalized to a peak value
%    of 1. It is called unitDensity for historical reasons.  The
%    normalization to peak of 1 is just a convention, the quantity that
%    matters is the product obj.density*ojb.unitDensity.  Note that in
%    other similar routines (e.g., Lens), we do not normalize the
%    obj.unitDensity to a peak of 1.
%
%    obj.density is the peak optical density (sometimes just called optical
%    density). The interpretation as peak optical density depends on the
%    convention followed here of normlizing obj.unitDensity to a peak of 1.
%    One estimate of the average (across observers) peak density ise 0.28,
%    with a range of 0.17 to 0.48. Our default is 0.35, matching that of
%    the underlying data file from CVRL.
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
%    Macular pigment peak optical density varies across the retina. Use the
%    eccDensity method to obtain peak optical density across the human
%    retina.
%
% Inputs:
%    None.
%
% Outputs:
%    obj - The created macular pigment object
%
% Optional key/value pairs:
%    name        - String. The object name. Default 'human macular'.
%    wave        - Vector. The wavelength support of returned functions.
%                  Default (400:10:700).
%    density     - Numeric. The peak optical density. Default 0.35. 
%    unitDensity - Vector. Underlying stored absorbance.  Should be normalized
%                  to a maximum of 1. Default is read in from ISETBio file
%                  macularPigment.mat, and normalized to max of 1. This
%                  combined with the default value for density leads to
%                  standard estimates. If this is passed, it should be on
%                  the same wavelength support as the wave (not wave_)
%                  property. If you set this when object is instantiated,
%                  the wave_ property is set to wave.
%

% History:
%    xx/xx/16  HJ   ISETBIO TEAM, 2016
%    02/15/18  jnm  Formatting
%    12/16/20  dhb  Mostly comments, but also changed name of set to allows
%                   setting of unitDensity_, not unitDensity. 
%    12/18/20  dhb  More comments.
%    12/20/20  npc  Made eccDensity an object method.
% Examples:
%{
    macobj = Macular;
    figure; plot(macobj.wave,macobj.transmittance);
    peripheralDensity = macobj.eccDensity(5);
%}
%{
    % Get and plot macular absorbance
    macobj = Macular;
    figure; plot(macobj.wave,macobj.unitDensity);

    % Set macular pigment absorbance to all ones. 
    % This gets set on the underlying wave
    % support.
    macobj.unitDensity = ones(size(macobj.wave));
    figure; plot(macobj.wave,macobj.unitDensity);
%}
%{
    % Create the object
    macobj = Macular;
    
    % Get and plot underlying macular absorbance on its wavelenght support.
    % This is splined onto the user's wavelength support when getting the
    % various properties that this object supports.
    figure; hold on; plot(macobj.wave_,macobj.unitDensity_,'r','LineWidth',1);

    % Get and plot macular absorbance as used on the 
    % users wavelength support, and for deriving other
    % properties.
    plot(macobj.wave,macobj.unitDensity,'ko','MarkerFaceColor','k','MarkerSize',4);

    % Set another wavelength support
    macobj.wave = 505:10:805;
    plot(macobj.wave,macobj.unitDensity,'go','MarkerFaceColor','g','MarkerSize',4);
%}

properties  % public properties
    % name - Name of this particular macular pigment object
    name;
    
    %density - macular pigment peak optical density
    density;
end

properties (SetObservable, AbortSet)
    %wave - wavelength samples in nm
    wave;
end

properties (SetAccess=private)
    %wave_ - The internal wavelength samples
    wave_;
    
    %unitDensity_ - unit density absorbance sampled with wave_ (this is
    %               usually just called absorbance). These values are
    %               interpolated to wave when a get is done on unitDensity.
    unitDensity_;
end

properties (Dependent)
    %unitDensity - spectral absorbance with unit pigment density (this is
    %              usually called just absorbance).
    unitDensity;

    %spectralDensity - unitDensity scaled by obj.density (this quantity
    %                  doesn't typically have its own name.
    spectralDensity;

    %transmittance - proportion of quanta transmitted
    transmittance;

    %absorptance - proportion of quanta absorbed
    absorptance;
end

methods  % public methods
    % constructor
    function obj = Macular(varargin)    
        p = inputParser;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('density', 0.35, @isscalar);
        p.addParameter('unitDensity', [], @isnumeric);
        p.addParameter('name', 'human lens', @isstr);
        p.parse(varargin{:});

        % set properties
        obj.name = p.Results.name;
        obj.wave = p.Results.wave(:);
        obj.density = p.Results.density;

        if isempty(p.Results.unitDensity)
            % Read and store default macular absorbance from file at its native
            % wavelength spacing of 390 to 830 nm.  This is splined onto the
            % wavelength support in the wave property upon get.
            %
            % The magic number 0.3521 below the maximum value in the file
            % being read, so dividing it out gives normalized ("unit")
            % density in the private variable we compute from.  This gets
            % put back when we compute "spectral" density, as the default
            % value for that is 0.35.  (Well, we lose the 0.0021 through
            % all this, but OK.)
            obj.wave_ = (390:830)';
            temp = ieReadSpectra('macularPigment.mat', ...
                obj.wave_);
            
            % You might think we would just normalize to the max of temp,
            % which is close to 0.3521.  But doing that will break our
            % validations because 0.3521 isn't exactly the max.  This then
            % produces small numerical differences in the validation checks
            % that we'd have to address, and that is a pain.
            obj.unitDensity_ = temp/0.3521;
        else
            % Store the passed unit density.  We assume it is on the
            % wavelength spacing in the wave_ property, and store it on
            % this spacing for splining out on read.
            obj.wave_ = obj.wave;
            obj.unitDensity_ = p.Results.unitDensity;
        end
    end

    % get methods for dependent variables
    function val = get.unitDensity(obj)
        % interpolate for wavelength samples
        %
        % Syntax:
        %   val = get.unitDensity(obj);
        %
        % Description:
        %    Retrieve the macular pigment object's unit density (aka
        %    absorbance)
        %
        % Inputs:
        %    obj - The macular pigment object
        %
        % Outputs:
        %    val - The unit density value
        %
        % Optional key/value pairs:
        %    None.
        
        val = interp1(obj.wave_, obj.unitDensity_, obj.wave, 'pchip');
        val = max(val, 0);
    end

    function val = get.spectralDensity(obj)
        % compute scaled absorbance
        %
        % Syntax:
        %   val = get.spectralDensity(obj);
        %
        % Description:
        %    Retrieve the macular pigments object's spectral density.
        %    This is the absorbance scaled by the peak optical density.
        %
        %    This is not standard terminology, as far as we know, and this
        %    quanity less likely to be useful than the absorptance or
        %    transmittance.
        %
        % Inputs:
        %    obj - The macular pigment object
        %
        % Outputs:
        %    val - The value of the object's spectral density
        %
        % Optional key/value pairs:
        %    None.

        val = obj.unitDensity * obj.density;
    end

    function val = get.transmittance(obj)
        % compute proportion of quanta transmitted
        %
        % Syntax:
        %   val = get.transmittance(obj);
        %
        % Description:
        %    Get the value of the macular pigment object's transmittance.
        %
        % Inputs:
        %    obj - The macular pigment object.
        %
        % Outputs:
        %    val - The transmittance value.
        %
        % Optional key/value pairs:
        %    None.
       
        val = 10 .^ (-obj.spectralDensity);
    end

    function val = get.absorptance(obj)
        % compute proportion of quanta absorbed
        %
        % Syntax:
        %   val = get.absorptance(obj);
        %
        % Description:
        %    Get the value of the macular pigment object's absorptance.
        %
        % Inputs:
        %    obj - The macular pigment object.
        %
        % Outputs:
        %    val - The absorptance value.
        %
        % Optional key/value pairs:
        %    None.

        val = 1 - obj.transmittance;
    end
    
    % set methods for dependent variables
    function set.unitDensity(obj, val)
        % Set underlying macular pigment absorbance.
        %
        % Syntax:
        %   set.unitDensity(obj, val)
        %
        % Description:
        %    Set the underlying macular pigment absorbance.  This should be
        %    on the wavelength spacing specified in property wave, not
        %    wave_.  But it is splined to obj.wave_ when it is stored.
        %
        % Inputs:
        %    obj - The Macular object
        %    val - The value to assign to unit density
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    None.

        if (length(val) ~= length(obj.wave))
            error('Must set absorbance (aka unitDensity) on wavelengths support of obj.wave_');
        end
        if (abs(max(val)-1) > 1e-2)
            % Can only check to precision of 1e-2 because of the values in
            % our historical validation files.
            error('Absorbance (aka unitDensity) should be normalized to a maximum of 1');
        end
        if (min(val) < -1e-4)
            error('Absorbance (aka unitDensity) should be be all positive');
        end
        
        % This would be better clipping to 0,1 after tighter checks above,
        % but not doing so to avoid breaking validtion files.
        val = max(val,0);
        obj.unitDensity_ = interp1(...
                obj.wave, val, obj.wave_, 'pchip');
    end
    
    % Method for computing macular pigment density as a function of
    % eccentricity
    density = eccDensity(obj, eccDeg, varargin);
end
    


end
