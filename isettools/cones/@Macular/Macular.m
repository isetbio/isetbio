classdef Macular < hiddenHandle
% Class for macular pigment properties
%
% Syntax:
%	macular = Macular('wave', wave)
%
% Description:
%    The human retina contains a pigment that covers the central (macular)
%    region. This macular pigment passes certain wavelengths of light more
%    than others. The pigment varies in density from central vision, where
%    it is highest, to increasingly peripheral vision.
%
%    This function returns several measures of the macular pigment
%    wavelength properties as a function of macular pigment density (high
%    in the fovea, lower in the near fovea).
%
%    The returned class object includes a variety of derived terms.
%    This should help to keep the relationship between entities straight.
%
%    Density is the estimated (average) peak density of the pigment across
%    a variety of observers.  They estimate the average (across observers)
%    peak density to be 0.28, with a range of 0.17 to 0.48.
%
%    Useful formulae
%       Absorbance spectra are normalized to a peak value of 1.
%       Absorptance spectra are the proportion of quanta actually absorbed.
%       Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% Inputs:
%    None required.
%
% Outputs:
%    The created macular pigment object
%
% Optional key/value pairs:
%    None.
%
% References:
%    * The original macular densities values were taken from the Stockman
%      site. Go to http://cvision.ucsd.edu, then click on Prereceptoral
%      filters.  At this point in time, I think the Psychtoolbox and the
%      new Stockman site are authoritative.

% History:
%    xx/xx/16  HJ   ISETBIO TEAM, 2016
%    02/15/18  jnm  Formatting

% Examples:
%{
    macular = Macular();
%}

properties  % public properties
    %density  macular pigment density
    density;
end

properties (SetObservable, AbortSet)
    %wave  wavelength samples in nm
    wave;
end

properties (Access=private)
    %wave_  The internal wavelength samples
    wave_;

    %unitDensity_  unit density absorbance sampled with wave_
    unitDensity_;
end

properties (Dependent)
    %unitDensity  spectral absorbance with unit pigment density
    unitDensity;

    %spectralDesnity  unitDensity scaled by obj.density
    spectralDensity;

    %transmittance  proportion of quanta transmitted
    transmittance;

    %absorptance  proportion of quanta absorbed
    absorptance;
end

methods  % public methods
    % constructor
    function obj = Macular(varargin)
        % Initialize defaults for Macular parameters
        %
        % Syntax:
        %   obj = Macular([varargin]);
        %
        % Description:
        %    Initialize the default values for the public properties: wave
        %    (400:10:700), density (.35), and unitDensity ([]). And then
        %    for the dependent and private object properties.
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
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('density', 0.35, @isscalar);
        p.addParameter('unitDensity', [], @isnumeric);

        p.parse(varargin{:});

        % set properties
        obj.wave = p.Results.wave(:);
        obj.wave_ = (390:830)';
        obj.density = p.Results.density;

        if isempty(p.Results.unitDensity)
            obj.unitDensity_ = ieReadSpectra('macularPigment.mat', ...
                obj.wave_) / 0.3521;
        else
            obj.unitDensity = p.Results.unitDensity;
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
        %    Retrieve the macular pigment object's unit density
        %
        % Inputs:
        %    obj - The macular pigment object
        %
        % Outputs:
        %    val - The unit Density value
        %
        % Optional key/value pairs:
        %    None.
        %
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
        %
        % Inputs:
        %    obj - The macular pigment object
        %
        % Outputs:
        %    val - The value of the object's spectral density
        %
        % Optional key/value pairs:
        %    None.
        %
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
        %
        val = 10 .^ (-obj.spectralDensity);
    end

    function val = get.absorptance(obj)
        % comptue proportion of quanta absorbed
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
        %
        val = 1 - obj.transmittance;
    end


    % set methods for dependent variables
    function set.unitDensity(obj, val)
        % interpolate for wavelength samples
        %
        % Syntax:
        %   set.unitDensity(obj, val)
        %
        % Description:
        %    Interpolate for the wavelength samples
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
        %
        obj.unitDensity_ = interp1(obj.wave, val,obj.wave_, 'pchip');
        obj.unitDensity_ = max(obj.unitDensity_, 0);
    end
end

methods (Static)
    function density = eccDensity(eccDeg)
        % Compute macular pigment optical density as funct. of eccentricity
        %
        % Syntax:
        %	density = macularDensity(eccDeg);
        %
        % Description:
        %    Compute the macular pigment's optical density as a function of
        %    its eccentricity.
        %
        % Inputs:
        %    eccDeg   - eccentricity in degrees
        %
        % Outputs:
        %	 density  - macular pigment optical density
        %
        % Notes:
        %   1) Macular pigment density is roughly symmetric and thus we
        %      approximate the 2D position by 1D eccentricity
        %   2) The lorentzian function is fitted from data grabbed from
        %      figure 2(B) in the reference paper. The data is stored
        %      in macularDensity.mat
        %
        % References:
        %   Putnam, C. M., & Bland, P. J. (2014). Macular pigment
        %   optical density spatial distribution measured in a subject
        %   with oculocutaneous albinism. Journal of Optometry, 7(4),
        %   241-245.
        %

        % Compute density with the lorentz function 
        % 
        % Here, we force the model to be symmetric and have 0 density
        % at infinite eccentricity
        density = 0.35 * 3.6028 ./ (eccDeg .^ 2 + 3.6028);
    end
end
end