function m = macularSet(m,param,val,varargin)
% Set parameters of the macular pigment structure
%
%     m = macularSet(m,param,val,varargin)
%
% See macularCreate for notes about the macular pigment in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Parameters
%
%   name
%   wave
%   unitDensity
%   density
%
% Examples:
%   m = macularCreate;
%   t = macularGet(m,'transmittance');
%   w = macularGet(m,'wave');
%   vcNewGraphWin; plot(w,t)
%
% Copyright ImagEval Consultants, LLC, 2013.

if notDefined('m'), error('Macular structure required'); end
if notDefined('param'), error('param required'); end
if notDefined('val'), error('val required'); end

%
param = ieParamFormat(param);

switch param
    case 'name'
        m.name = val;
    case {'wave', 'wavelengths', 'wavelength'}
        % change the sampling wavelength
        % this requires re-sampling of the data
        if ~isempty(varargin)
            % passed in extraval in varargin
            extrapval = varargin{1};
        else
            extrapval = 0;
        end
        
        val = val(:); % Make sure it's column vector
        wave = macularGet(m, 'wave');
        unitDensity = macularGet(m, 'unit density');
        m.unitDensity = interp1(wave, unitDensity, val, [], extrapval);
        m.wave = val;
        
    case 'unitdensity'
        % Spectral density
        % Accept wavelength in varargin{1}
        if isempty(varargin)
            % Checking for lengths
            assert(length(val) == length(macularGet(m, 'wave')), ...
                'Val should have same length as macular wavelength');
            m.unitDensity = val;
        else
            wave = varargin{1};
            wave = wave(:); % Make sure a column vector
            assert(length(val)==length(wave),'val and wave size mismatch');
            m.unitDensity = interp1(wave, val, m.wave, [], 0);
        end
    case 'density'
        % Density for this case
        m.density = val;
    otherwise
        error('Unknown parameter %s\n',param);
end

end