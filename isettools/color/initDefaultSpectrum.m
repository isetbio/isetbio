function ieStruct = initDefaultSpectrum(ieStruct, spectralType, wave)
% Create a wavelength spectrum structure and attach it to an ISET object
%
% Syntax:
%   ieStruct = initDefaultSpectrum(ieStruct, spectralType, [wave])
%
% Description:
%    Initialize a spectrum structure field of a structure. The spectrum
%    structure specifies the the  spectrum and sample wavelengths.
%
%    This is used for early code where "objects" were in fact structures
%    (e.g. scene, oi), and is primarly used internally to initialize
%    spectral fields.
%
%    The initialization consists only of attaching a field to the passed
%    structure that is called 'spectrum', and adding the 'wave' field to
%    that structure. Something else is required to add in an actual field
%    that specifies a spectrum.
%
%    This file contains examples of usage. To acces, type 'edit
%    initDefaultSpectrum.m' into the Command Window.
%
% Inputs:
%    ieStruct     - Struct. The ISET structure to which you wish to attach
%                   the created spectrum structure.
%    spectralType - String. One of the (currently) supported types below.
%          multispectral: Default 400:10:700 nm
%          monochrome: 550 nm
%          custom: The user supplies the wavelength samples (deprecated)
%    wave         - (Optional) Vector. The wavelength samples required for
%                   a custom spectral type (should be deprecated)
%
% Outputs:
%    ieStruct     - Struct. The structure with the spectrum field added.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [NOTE - DHB: This was very early iset code, and I'm not really sure
%      how it is used or whether we really want it going forward. But it
%      is called in many places, so it won't be easy to get rid of.]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    10/27/17  jnm  Comments & formatting
%    07/10/19  JNM  Formatting update

% Examples:
%{
    scene = sceneCreate;
    scene = initDefaultSpectrum(scene, 'monochrome');
    scene = initDefaultSpectrum(scene, 'multispectral');
    scene = initDefaultSpectrum(scene, 'custom', 400:50:700);
%}

if notDefined('ieStruct'), error('Input ieStruct required.'); end
if notDefined('spectralType'), spectralType = 'hyperspectral'; end

switch lower(spectralType)
    case {'multispectral', 'hyperspectral','spectral'}
        ieStruct.spectrum.wave = (400:10:700)';

    case 'monochrome'
        ieStruct.spectrum.wave = 550;

    case 'custom'
        if notDefined('wave')
            error('wave required for custom spectrum');
        end
        ieStruct.spectrum.wave = wave(:);

    otherwise
        error('spectralType not yet defined.');
end

end
