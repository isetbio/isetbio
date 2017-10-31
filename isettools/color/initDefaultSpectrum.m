function object = initDefaultSpectrum(object, spectralType, wave)
% Create a wavelength spectrum structure and attach it to an ISET object
%
% Syntax:
%   object = initDefaultSpectrum(object, spectralType, wave)
%
% Description:
%    The spectrum structure specifies the sample wavelengths.
%
%    We use only three spectral types at present.  These are
% Inputs:
%    spectralType - One of the (currently) three supported types below.
%           Multispectral - 400:10:700 nm
%           Monochrome    - 550 nm
%           Custom        - The user supplies the wavelength samples
%    wave         - (Optional) the wavelength samples required for a custom
%                   spectral type
%
% Outputs:
%    obj          - the ISET object to which you wish to attach the created
%                   spectrum structure
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    10/27/17  jnm  Comments & formatting

% Examples:
%{
   scene = sceneCreate;
   scene = initDefaultSpectrum(scene, 'monochrome');
   scene = initDefaultSpectrum(scene, 'multispectral');
   scene = initDefaultSpectrum(scene, 'custom', 400:50:700);
%}

if notDefined('object'), error('Object required.'); end
if notDefined('spectralType'), spectralType = 'hyperspectral'; end

switch lower(spectralType)
    case {'spectral', 'multispectral', 'hyperspectral'}
        object.spectrum.wave = (400:10:700)';
        
    case 'monochrome'
        object.spectrum.wave = 550;
        
    case 'custom'
        if notDefined('wave')
            error('wave required for custom spectrum');
        end
        object.spectrum.wave = wave(:);
        
    otherwise
        error('spectralType not yet defined.');
end

end
