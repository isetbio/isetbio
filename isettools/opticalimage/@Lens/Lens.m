classdef Lens < handle
% Class for human lens pigment properties
%
% Syntax:
%   lens = LENS()
%
% Description:
%    returns a lens object whose public properties are the lens name,
%    wave, and density.  The density refers to the density of the lens
%    pigment.  
%
%    The lens object stores the calibrated lens density values in
%    private variables (obj.wave_ and obj.unitDensity_). These values
%    were taken from the PsychToolbox and/or the Stockman site. Go to
%    http://cvision.ucsd.edu, then click on Prereceptoral filters.
%
%    Other values, such as absorbance, absorptance and transmittance,
%    are derived from these parameters.
%
%    Useful factoids:
%
%     * Absorbance spectra are normalized to a peak value of 1.
%     * Absorptance spectra are the proportion of quanta actually absorbed.
%     
%    Equations: 
%         absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%         transmittance      = 1 - absorptance 
%
%    There are examples contained in the code. To access, type 'edit
%    lens.m' into the Command Window.
%
% Inputs:
%    None required.
%
% Outputs:
%    The lens object
%
% Optional key/value pairs:
%    **Needs to be filled out**

% History:
%    xx/xx/13  HJ/BW  ISETBIO Team 2013.
%    03/05/18  jnm    Formatting

% Examples:
%{
    lens = Lens();
%}
%{
    thisLens = Lens('wave', 400:10:700, 'density', 1, 'name', 'my lens');
%}
%{
    lens = Lens();
    t = lens.transmittance;
    a = lens.absorptance;
    vcNewGraphWin; 
    plot(lens.wave,t,'r-',lens.wave,a,'g-',lens.wave,a+t,'k--');
    xlabel('Wave (nm)'); ylabel('Fraction');
    legend({'transmittance','absorptance','sum'});
%}

properties  % public properties
    % name - Name of this particular lens object
    name;

    % density - pigment density
    density;

    % wave - wavelength samples in nm
    wave;
end

properties (Access=private)
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
        %
        % Syntax:
        %   thisLens = 
        %       Lens('wave', 400:10:700, 'density', 1, 'name', 'my lens');
        %

        % parse input
        p = inputParser;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('density', 1, @isscalar);
        p.addParameter('unitDensity', [], @isnumeric);
        p.addParameter('name', 'human lens', @isstr);

        p.parse(varargin{:});

        obj.name = p.Results.name;

        % set properties
        obj.wave = p.Results.wave(:);
        obj.wave_ = (390:830)';
        obj.density = p.Results.density;

        if isempty(p.Results.unitDensity)
            obj.unitDensity_ = ieReadSpectra('lensDensity.mat', ...
                obj.wave_);
        else
            obj.unitDensity = p.Results.unitDensity;
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
        %    obj - The lens object
        %
        % Outputs:
        %    val - The unit density
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
        %    obj - The lens object
        %
        % Outputs:
        %    val - The spectral density
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
        %    obj - The lens object
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
        %         name                      - String. The Lens object name
        %         {wave, wavelength}        - Numerical. Wavelength vector
        %         {absorbance, unitDensity} - Numerical. Length match wave
        %         density                   - Numerical. Scalar value.
        %    val   - The value to assign
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
                assert(length(val) == length(lens.wave), ...
                    'Val should have same length as lens wavelength');
                obj.unitDensity_ = interp1(...
                    obj.wave, val, obj.wave_, 'pchip');
                obj.unitDensity_ = max(obj.unitDensity, 0);
            case 'density'
                assert(isscalar(val), 'val should be scalar');
                obj.density = val;
            otherwise
                error('Unknown parameter %s\n', param);
        end
    end

    % set methods for dependent variables
    %         function obj = set.unitDensity(obj, val)
    %             % interpolate for wavelength samples
    %             obj.unitDensity_ = interp1(...
    %                 obj.wave, val, obj.wave_, 'pchip');
    %             obj.unitDensity_ = max(obj.unitDensity, 0);
    %         end
    %
    %         function set.wave(obj, val)
    %             obj.wave = val;
    %         end
    %
    %         function set.density(obj, val)
    %             obj.density = val;
    %         end

end
end