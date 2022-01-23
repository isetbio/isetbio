classdef receptorPigment < hiddenHandle
% Abstract class (interface) for photopigment objects:
%
% This class contains properties for the photopigment absorption properties
% of a single cone cell. For the full cone mosaic, see the cMosaic class.
%
% Subclasses of @receptorPigment inherit spectral properties and methods
% from @receptorPigment.  Using this superclass scheme, enforces
% consistency amongst all subclasses in terms of how they handle spectral
% properties.
% 
% The subclasses define their own aperture geometrical properties and
% methods. For example, the old  @photoPigment has rectangular-shaped
% aperture whereas the new @cPhotopigment has disk-shaped aperture.
%
% Description:
%    Our understanding of the terminology and best conventions for
%    describing this type of data has evolved over time, and the property
%    names differ somewhat from the conventions we would adopt today.
%    Changing the names in the code will not be backwards compatibile, so
%    we have done our best to comment and explain here.
%
%    obj.absorbance  - the absorbance spectra of the L, M, and S cones,
%    each normalized to a peak value of 1. Values are in the columns, with
%    a separate column for L, M and S. The column-wise arrangement applies
%    to this and to the other L, M, and S spectra below. This quantity is
%    called obj.unitDensity in the Lens and Macular objects, and in the
%    Lens object it is not normalized. The normalization to peak of 1 is
%    just a convention, the quantity that matters is the product
%    obj.opticalDensity*ojb.absorbance, and there are those who might call
%    that product the absorbance.
%
%    obj.opticalDensity - the peak optical density (sometimes just called
%    optical density). The interpretation as peak optical density depends
%    on the convention followed here of normlizing obj.absorbance to a peak
%    of 1. There are three entries to this vector, one each for the L, M,
%    and S cones.
%
%    obj.peakEfficiency is the probability that a photopigment absorption
%    leads to an isomerization.  There are three entries to this vector,
%    one each for the L, M, and S cones.  This property is unfortunately
%    named, as is it isn't the peak of anything.
%
%    obj.absorptance is the absorptance spectrum.  This tells us the
%    probability that a photon of a given wavelength is absorbed as it
%    passes through a layer of photopigment with total absorbance given
%    by obj.opticalDensity*obj.absorbance.
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
%       
%           absorptance = 1 - 10.^(-opticalDensity * absorbance). 
%
%       Absorbance spectra are normalized to a peak value of 1, and then
%       scaled by optical density to get the not normalized absorbance.
%
%       Absorptance spectra are the proportion of quanta actually absorbed.
%       This is the term used in this routine.
%       
%       In this routine, again for historical reasons, opticalDensity is
%       just called density.  In the literature, this is sometimes called
%       peak optical density.       
%
%    The absorbance data that drive this routine are stored on wavelength
%    support in property wave_ in property absorbance.  Typically wave_ is
%    set to a large wavelength support and then interpolated onto the
%    support in propety wave.  You can set absorbance after the object is
%    instantiated, but you can't change wave_.  When you set absorbance, it
%    should be on wavelength support wave, and it is splined onto the
%    wavelength support in wave_ before being stored in obj.absorbance.
%
%    This design does not prevent you from setting absorbance on wavelength
%    support very different from that being previously used to store the
%    data, which could lead to extrapolation errors.  To avoid this, if you
%    want to use custom data, you may be better off creating the object
%    with the desired data on the wavelength support you intend to use.
%    That said, the default values are read in and stored on wave_ support
%    of 390:830 at 1 nm spacing, which is good for most applications.
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
%

% History:
%   12/19/2020  npc   Wrote it from @photoPigment
%

% Public properties 
properties 
    % opticalDensity - photopigment optical densities for different cone types
    opticalDensity;

    % peakEfficiency - peak absorptance efficiency
    peakEfficiency;
    
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
    
    % quantalEfficiency - probability of isomerization in quantal units.
    % Gives the probability that an incident photon isomerizes
    % photopigment. These are actually useful.  Does not take into account
    % effect of inert pigments (lens, macular pigment).
    quantalEfficiency;
end

% Protected properties. All receptorPigment subclasses can read these, but they cannot set them. 
properties (SetAccess = protected)
end

% ReceptorPigment subclasses do not have any access to these properties directly
properties (SetAccess = protected, GetAccess = protected)
        
end

properties (SetObservable, AbortSet)
    % wave - wavelength samples
    wave;
end

properties (SetAccess = public)
    % I made these public so I could change them from a script.  But NC may
    % want us to do this a different way with set operations.
    %
    
    % wave_ - The internal wavelength samples
    wave_;

    % absorbance_ - The absorbance data sampled at wave_
    absorbance_;
end

% Abstract, public methods. Each subclass of @receptorPigment *must* implenent
%  its own version of all functions listed as abstract. If it does not, 
% it cannot instantiate objects.
methods(Abstract)
    description(obj, varargin)
end

% Public methods
methods
    % Constructor
    function obj = receptorPigment(varargin)
        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
        p.addParameter('absorbance', [], @isnumeric);
        p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
        p.parse(varargin{:});
        
        % set object properties
        obj.wave = p.Results.wave(:); 
        obj.opticalDensity = p.Results.opticalDensity(:);
        obj.peakEfficiency = p.Results.peakEfficiency(:);

        % If absorbance is not specified, we obtain it using the defaults
        % of coneAbsorbanceReadData.
        % 
        if isempty(p.Results.absorbance)
            obj.wave_ = (390:830)';
            obj.absorbance_ = coneAbsorbanceReadData(p.Unmatched, ...
                'wave', obj.wave_);
        else
            obj.wave_ = p.Results.wave;
            obj.absorbance_ = p.Results.absorbance;
        end
    end % Constructor
    
    % Getter methods for dependent variables
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
    
end % Public methods


end