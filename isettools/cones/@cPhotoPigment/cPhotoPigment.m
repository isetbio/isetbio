classdef cPhotoPigment < hiddenHandle
% Class for single cone photopigment and related properties
%
% Syntax:
%	pigment = cPhotoPigment;
%
% Description:
%    This class contains properties for the photopigment absorption
%    properties of a single cone cell. This class is derived from the 
%    @photoPigment, which assumes a rectangular shaped cone, as does the
%    @coneMosaic class. @cPhotoPigment assumes a circular shaped cone and
%    is supposed to work with the @cMosaic class.
%
%    For the full cone mosaic, see the @cMosaic class
%
%    Most of the terms represented here are descriptions of the
%    photopigment itself. In addition, there are a few terms the
%    capture the effective optical size of the photopigment absorption.
%
%    Default parameters are determined by underlying routines that get
%    the required data types. Unmatched key/value pairs passed to
%    photoPigment are passed on to some underlying routines and can be
%    used to adjust the parameters obtained. See help for each routine
%    for what the available key/value pairs are.
%
% Input:
%	 None required.
%
% Output:
%    pigment          - The created cPhotoPigment object.
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
%    'diameter'       - Diameter of the light collecting aperture Default 2e-6.
%
% See Also:
%   photoPigment, cMosaic, Macular, Lens
%

% History:
%    12/07/20  NPC Wrote it by adapting photoPigment

properties  % public properties
    % opticalDensity - photopigment optical densities for cones
    opticalDensity;

    % peakEfficiency - peak absorptance efficiency
    peakEfficiency;

    % diameter of the light collecting cone aperture in meters
    diameter;
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

    % area - The area of the light collecting aperture.
    area;
end

properties(Access = private)  % private properties
    % wave_ - The internal wavelength samples
    wave_;

    % absorbance_ - The absorbance data sampled at wave_
    absorbance_;
end

methods  % public methods
    % constructor
    function obj = cPhotoPigment(varargin)

        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
        p.addParameter('absorbance', [], @isnumeric);
        p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
        p.addParameter('diameter', (4.0/sqrt(pi))*1e-6, @isnumeric);
        p.parse(varargin{:});

        % set object properties
        obj.wave = p.Results.wave(:);
        obj.wave_ = (390:830)';
        obj.opticalDensity = p.Results.opticalDensity(:);
        obj.peakEfficiency = p.Results.peakEfficiency(:);
        obj.diameter = p.Results.diameter;
        
        % Assert that property dimensions are consistent
        assert(numel(obj.opticalDensity) == numel(obj.peakEfficiency), ...
            sprintf('optical density dimensionality does not match that of peak efficiency'));

        % If absorbance is not specified, we obtain it using the defaults
        % of coneAbsorbanceReadData. 
        if isempty(p.Results.absorbance)
            obj.absorbance_ = coneAbsorbanceReadData(p.Unmatched, ...
                'wave', obj.wave_);
        else
            obj.absorbance = p.Results.absorbance;
        end
        
        % Assert that property dimensions are consistent
        if (~isempty(p.Results.absorbance))
            assert(numel(obj.opticalDensity) == size(obj.absorbance,2), ...
                sprintf('optical density dimensionality does not match that of absorbance'));  
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

    function val = get.absorptance(obj)
        % Compute the absorptances
        %
        % Syntax:
        %   obj = get.absorptance(obj)
        %
        % Description:
        %    Compute the absorptances
        %
        % Inputs:
        %    obj - The cPhotoPigment object
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
        % Compute the quanta fundamentals
        %
        % Syntax:
        %   obj = get.quantaFundamentals(obj)
        %
        % Description:
        %    Compute and return the quanta fundamentals for the photo
        %    pigment object obj.
        %
        % Inputs:
        %    obj - The cPhotoPigment object
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
        % Compute the energy fundamentals
        %
        % Syntax:
        %   obj = get.energyFundamentals(obj)
        %
        % Description:
        %    Retrieve the energy fundamentals from the cPhotoPigment object
        %
        % Inputs:
        %    obj - The cPhotoPigment object
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
        % Compute  the light-collecting area 
        %
        % Syntax:
        %   obj = get.area(obj)
        %
        % Description:
        %    Compute  the light-collecting area 
        %
        % Inputs:
        %    obj - The cPhotoPigment object
        %
        % Outputs:
        %    val - The area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = pi * (0.5*obj.diameter)^2;
    end


    % set method for dependent variable
    function set.absorbance(obj, val)
        % Set the the absorbance values
        %
        % Syntax:
        %   obj = set.absorbance(obj, val)
        %
        % Description:
        %    Set the absorbance value
        %
        % Inputs:
        %    obj - The cPhotoPigment object
        %    val - The absorbance value to set
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    None.
        %
        
        obj.absorbance_ = interp1(obj.wave, val, obj.wave_, ...
            'linear', 'extrap');
        obj.absorbance_ = ieClip(obj.absorbance_, 0, 1);
    end
end

methods (Static)
    % When we have them, they go here
end
end