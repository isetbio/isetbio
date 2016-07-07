function val = lensGet(lens,param,varargin)
% Get parameters of the lens structure
%
%     val = lensGet(lens,param,varargin)
%
% See lensCreate for notes about the properties in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Parameters
%
%   name
%   type          - 'lens'
%   unitDensity   - Read in from lensDensity.mat file, based on Sharp
%   density       - single value
%   transmittance
%   absorbance
%   absorptance (absorption)
%
%   Absorbance spectra are normalized to a peak value of 1.
%   Absorptance spectra are the proportion of quanta actually absorbed.
%   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% Examples:
%  lens = lensCreate; w = lensGet(lens,'wave');
%  vcNewGraphWin; plot(w,lensGet(lens,'absorptance'))
%  hold on; plot(w,lensGet(lens,'transmittance'))
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('lens'), error('Lens structure required'); end
if notDefined('param'), error('param required'); end

%
param = ieParamFormat(param);

switch param
    case 'name'
        val = lens.name;
    case 'type'
        val = lens.type;
    case 'wave'
        val = lens.wave;
        
        
    case {'absorbance','unitdensity'}
        % We accept user defined sampling wavelength in varargin{1}
        if isempty(varargin)
            val = lens.unitDensity;
        else
            wave = lensGet(lens, 'wave');
            val = interp1(wave, lens.unitDensity, varargin{1});
        end
    case 'density'
        % Assumed density for this instance
        val = lens.density;

    case {'spectraldensity'}
        % Unit density times the density for this structure
        if isempty(varargin)
            u = lensGet(lens,'unit density');
        else
            u = lensGet(lens,'unit density', varargin{1});
        end
        d = lensGet(lens,'density');
        val = u*d;
        
    case 'transmittance'
        % Proportion of quanta transmitted
        if isempty(varargin)
            val = 10.^(-lensGet(lens,'spectral density'));
        else
            val = 10.^(-lensGet(lens,'spectral density', varargin{1}));
        end
        
    case {'absorptance','absorption'}
        % Proportion of quanta absorbed
        if isempty(varargin)
            val = 1 - 10.^(-lensGet(lens,'spectral density'));
        else
            val = 1 - 10.^(-lensGet(lens,'spectral density', varargin{1}));
        end
        
    otherwise
        error('Unknown parameter %s\n',param);
end

return