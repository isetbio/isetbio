function val = get(lens,param,varargin)
% Get parameters of the lens structure
%
%     val = get(lens,param,varargin)
%
% See Lens object for notes about the properties in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Parameters
%   name                     - this lens name
%   {absorbance,unitDensity} - from lensDensity.mat file, based on Sharp
%   density                  - lens pigment density
%   transmittance            - 10^(-(spectral density))
%   absorptance (absorption) - 1 - transmittance;
%
% Note:
%  Absorbance spectra are normalized to a peak value of 1.
%  Absorptance spectra are the proportion of quanta actually absorbed.
%  Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% HJ/BW, ISETBIO TEAM, 2016

% parse parameters
p = inputParser;
p.addRequired('param',@isstr);
p.addParameter('wave',lens.wave,@isvector);

p.parse(param,varargin{:});
wave = p.Results.wave;

% set wavelength
lens = lens.copy;
lens.wave = wave;

switch ieParamFormat(param)
    case 'name'
        val = lens.name;
    case 'wave'
        val = lens.wave;
    case {'absorbance','unitdensity'}
        val = lens.unitDensity;
    case 'density'
        val = lens.density;
    case 'spectraldensity'
        % Unit density times the density for this structure
        val = lens.spectralDensity;
    case 'transmittance'
        % Proportion of quanta transmitted
        val = lens.transmittance;
    case {'absorptance','absorption'}
        % Proportion of quanta absorbed
        val = lens.absorptance;
    otherwise
        error('Unknown parameter %s\n',param);
end

end