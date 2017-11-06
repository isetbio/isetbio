function ieStruct = initDefaultSpectrum(ieStruct, spectralType, wave)
% Create a wavelength spectrum structure and attach it to an ISET object
%
% Syntax:
%   ieStruct = initDefaultSpectrum(ieStruct, spectralType, [wave])
%
% Description:
%    Initialize a spectrum structure field of a structure.  The spectrum
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
% Inputs:
%    ieStruct     - The ISET object to which you wish to attach the created
%                   spectrum structure
%    spectralType - One of the (currently) three supported types below.
%           'multispectral' - 400:10:700 nm. The strings 'spectral' and
%                             'hyperspectral' are synonyms for this.
%           'monochrome'    - 550 nm
%           'custom'        - The user supplies the wavelength samples
%    wave         - (Optional) the wavelength samples required for a custom
%                   spectral type
%
% Outputs:
%    ieStruct     - The structure with the spectrum field added.
%
% Notes:
%   * [NOTE - DHB: This was very early iset code, and I'm not really sure
%      how it is used or whether we really want it going forward.  But it
%      is called in many places, so it won't be easy to get rid of.]
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

if notDefined('ieStruct'), error('Input ieStruct required.'); end
if notDefined('spectralType'), spectralType = 'hyperspectral'; end

switch lower(spectralType)
    case {'spectral', 'multispectral', 'hyperspectral'}
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
