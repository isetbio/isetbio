function lens = lensSet(lens,param,val,varargin)
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
%
%   Absorbance spectra are normalized to a peak value of 1.
%   Absorptance spectra are the proportion of quanta actually absorbed.
%   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% Examples:
%  lens = lensCreate; w = lensGet(lens,'wave');
%  vcNewGraphWin; plot(w,lensGet(lens,'absorptance'))
%  hold on; plot(w,lensGet(lens,'transmittance'))
%  lens = lensSet(lens,'density',0.5);
%  vcNewGraphWin; plot(w,lensGet(lens,'absorptance'))
%  hold on; plot(w,lensGet(lens,'transmittance'))
%
% Copyright ImagEval Consultants, LLC, 2005.

if      notDefined('lens'), error('Lens structure required');
elseif ~isequal(lens.type,'lens'), error('Not a lens structure');
end

if notDefined('param'), error('param required'); end
if notDefined('val'),   error('val required');   end


param = ieParamFormat(param);

switch param
    case 'name'
        lens.name = val;
    case {'wave', 'wavelength'}
        % Sampling wavelength
        % This will require the interpolation on the data
        if ~isempty(varargin)
            % passed in extraval in varargin
            extrapval = varargin{1};
        else
            extrapval = 0;
        end
        wave = lensGet(lens, 'wave');
        unitDensity = lensGet(lens, 'unit density');
        
        val = val(:); % Make sure it's a column vector
        lens.unitDensity = interp1(wave, unitDensity, val, [],extrapval);
        lens.wave = val;
    case {'absorbance','unitdensity'}
        % This might not be the unit density, which has me bummed.  We
        % should deal with this in lensCreate.
        %
        % ieReadSpectra('lensDensity.mat',wave);
        % See macularCreate.
         % Spectral density
        % Accept wavelength in varargin{1}
        if isempty(varargin)
            % Checking for lengths
            assert(length(val) == length(lensGet(lens, 'wave')), ...
                'Val should have same length as lens wavelength');
            lens.unitDensity = val;
        else
            wave = varargin{1};
            wave = wave(:); % Make sure a column vector
            assert(length(val)==length(wave),'val and wave size mismatch');
            lens.unitDensity = interp1(wave, val, lens.wave, [], 0);
        end
    case 'density'
        % Assumed density for this instance
        assert(isscalar(val), 'val should be scalar');
        lens.density = val;      
    otherwise
        error('Unknown parameter %s\n',param);
end

end