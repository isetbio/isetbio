function val = get(lens,param,varargin)
% Get parameters of the lens structure
%
%     val = get(lens,param,varargin)
%
% See lensCreate for notes about the properties in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Parameters
%
%   name                     - this lens name
%   {absorbance,unitDensity} - Read in from lensDensity.mat file, based on Sharp
%   density                  - single value
%   transmittance            - 10^(-(spectral density))
%   absorptance (absorption) - 1 - transmittance;
%
%  Absorbance spectra are normalized to a peak value of 1.
%  Absorptance spectra are the proportion of quanta actually absorbed.
%  Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% Examples:
%  lens = lensCreate; w = lensGet(lens,'wave');
%  vcNewGraphWin; plot(w,lensGet(lens,'absorptance'))
%  hold on; plot(w,lensGet(lens,'transmittance'))
%
% Copyright ImagEval Consultants, LLC, 2005.

p = inputParser;
p.addRequired('param',@isstr);
p.addParameter('wave',lens.wave,@isvector);

p.parse(param,varargin{:});
param = p.Results.param;
wave = p.Results.wave;

lens = lens.copy;
lens.wave = wave;

%
param = ieParamFormat(param);

switch param
    case 'name'
        val = lens.name;
    case 'wave'
        val = lens.wave;
        
    case {'absorbance','unitdensity'}
        % We accept user defined sampling wavelength in varargin{1}
        val = lens.unitDensity;
    case 'density'
        % Assumed density for this instance
        val = lens.density;

    case {'spectraldensity'}
        % Unit density times the density for this structure
        u = lens.get('unit density');
        d = lens.get('density');
        val = u*d;
        
    case 'transmittance'
        % Proportion of quanta transmitted
        val = 10.^(-lens.get('spectral density'));

    case {'absorptance','absorption'}
        % Proportion of quanta absorbed
        val = 1 - 10.^(-lens.get('spectral density'));
        
    otherwise
        error('Unknown parameter %s\n',param);
end

end
