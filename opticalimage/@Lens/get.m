function val = get(lens, param, varargin)
% Get parameters of the lens structure
%
% Syntax:
%   val = get(lens, param, [varargin])
%
% Description:
%    See Lens object for notes about the properties in general and the
%    formulae relating absorbance and absorptance and transmittance.
%
%    * Absorbance spectra are normalized to a peak value of 1.
%    * Absorptance spectra are the proportion of quanta actually absorbed.
%    * Equation: absorptanceSpectra = 1 - 10 .^ (-OD * absorbanceSpectra)
%
% Inputs:
%    lens     - Object. The lens object to retrieve information from
%    param    - String. The parameter which you want to know about. Options
%               for the param, and the associated return type include:
%       name            - String. The lens name.
%       wave            - Vector. Wavelength samples in nm.
%       {absorbance, unitDensity}
%                       - File. The unitDensity information. Default is the
%                         file lensDensity.mat based on Sharpe and Stockman
%       density         - Numeric. Lens pigment density.
%       transmittance   - Numeric. The transmittance, calculated by
%                         10 ^ (-(spectral density)).
%       {absorptance, absorption}
%                       - Numeric. The absorptance value, calculated by
%                         1 - transmittance;
%       spectraldensity - Numeric. Unit density times the density.
%
% Outputs:
%    val      - VARIES. The value of the requested parameter. See param's
%               listed options above for type information.
%
% Optional key/value pairs:
%    None.
%
% See also:  opticsSet, opticsGet

% History:
%    xx/xx/16  HJ/BW  ISETBIO TEAM, 2016
%    03/05/18  jnm    Formatting
%    07/03/19  JNM    Formatting update

% parse parameters
p = inputParser;
p.addRequired('param', @isstr);
p.addParameter('wave', lens.wave, @isvector);

p.parse(param, varargin{:});
wave = p.Results.wave;

% set wavelength
lens.wave = wave;

switch ieParamFormat(param)
    case 'name'
        val = lens.name;
    case 'wave'
        val = lens.wave;
    case {'absorbance', 'unitdensity'}
        val = lens.unitDensity;
    case 'density'
        val = lens.density;
    case 'spectraldensity'
        % Unit density times the density for this structure
        val = lens.spectralDensity;
    case 'transmittance'
        % Proportion of quanta transmitted
        val = lens.transmittance;
    case {'absorptance', 'absorption'}
        % Proportion of quanta absorbed
        val = lens.absorptance;
    otherwise
        error('Unknown parameter %s\n', param);
end

end