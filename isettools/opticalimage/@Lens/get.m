function val = get(lens, param, varargin)
% Get parameters of the lens structure
%
% Syntax:
%   val = get(lens, param, varargin)
%
% Description:
%    See Lens object for notes about the properties in general and the
%    formulae relating absorbance and absorptance and transmittance.
%
% Inputs:
%    lens     - Object. The lens object to retrieve information from
%    param    - String. The parameter which you want to know about. Options
%               for the param include:
%       name                      - The lens name
%       {absorbance, unitDensity} - These will be retrieved from
%                                   lensDensity.mat file, based on Sharp.
%       density                   - Lens pigment density
%       transmittance             - The transmittance, calculated by
%                                   10 ^ (-(spectral density))
%       absorptance (absorption)  - The absorptance value, calculated by
%                                   1 - transmittance;
%    varargin - (Optional) Additional information that may be required in
%               the call, such as unit type.
%
% Outputs:
%    val      - The value of the requested parameter.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Absorbance spectra are normalized to a peak value of 1.
%    * Absorptance spectra are the proportion of quanta actually absorbed.
%    * Equation: absorptanceSpectra = 1 - 10 .^ (-OD * absorbanceSpectra)
%

% History:
%    xx/xx/16  HJ/BW  ISETBIO TEAM, 2016
%    03/05/18  jnm    Formatting

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