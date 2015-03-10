function val = macularGet(m,param,varargin)
% Get parameters of the macular pigment structure
%
%     val = macularGet(m,param,varargin)
%
% See macularCreate for notes about the macular pigment in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Parameters
%
%   name
%   type          - 'macular'
%   unitDensity   - Read in from macularPigment.mat file, based on Sharp
%   density       - single value
%   transmittance
%   absorbance
%   absorptance (absorption)
%
%
% Examples:
%   m = macularCreate; w= macularGet(m,'wave');
%   vcNewGraphWin; plot(w,macularGet(m,'absorbance'))
%   hold on; plot(w,macularGet(m,'transmittance'))
%
% Copyright ImagEval Consultants, LLC, 2013

if notDefined('m'), error('Macular structure required'); end
if notDefined('param'), error('param required'); end

%
param = ieParamFormat(param);

switch param
    case 'name'
        val = m.name;
    case 'type'
        val = m.type;
    case 'wave'
        val = m.wave;
    
    %   Absorbance spectra are normalized to a peak value of 1.
    %   Absorptance spectra are the proportion of quanta actually absorbed.
    %   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
    case {'absorbance', 'unitdensity'}
        % This is defined by Sharp, 1999.  To load use
        % ieReadSpectra('macularPigment.mat',wave);
        % See macularCreate.
        % We accept user defined wavelength in varargin
        if isempty(varargin)
            val = m.unitDensity;
        else
            wave = varargin{1};
            val = interp1(macularGet(m, 'wave'), m.unitDensity, wave);
        end
        
    case 'density'
        % Assumed density for this instance
        % This should be single scaler
        val = m.density;

    case {'spectraldensity'}
        % Unit density times the density for this structure
        if isempty(varargin)
            u = macularGet(m, 'unit density');
        else
            u = macularGet(m, 'unit density', varargin{1});
        end
        d = macularGet(m, 'density');
        val = u * d;
        
    case 'transmittance'
        % Proportion of quanta transmitted
        % accpet user defined wavelength in varargin{1}
        if isempty(varargin)
            val = 10.^(-macularGet(m, 'spectral density'));
        else
            val = 10.^(-macularGet(m, 'spectral density', varargin{1}));
        end
        
    case {'absorptance','absorption'}
        % Proportion of quanta absorbed
        if isempty(varargin)
            val = 1 - macularGet(m, 'transmittance');
        else
            val = 1 - macularGet(m, 'transmittance', varargin{1});
        end
        
    otherwise
        error('Unknown parameter %s\n',param);
end

return


