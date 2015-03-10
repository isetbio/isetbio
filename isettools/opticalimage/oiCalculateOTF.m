function [otf,fSupport] = oiCalculateOTF(oi,wave,unit)
%Calculate the optical transfer function for the optical image
%
%   [otf,fSupport] = oiCalculateOTF(oi,[wave],[unit = 'cyclesPerDegree'])
%
% The optical transfer function (OTF) is derived from the optics parameters
% of an optical image.  The frequency units are cycles per degree by
% default.  However, by setting the variable unit='millimeter' or 'micron'
% the frequency units can be changed to cycles/{millimeter,micron}. 
%
% This routine is used for diffraction limited, shift-invariant, or the
% human otf. Calculations of OTF for ray trace dat are handled in the ray
% trace routines, rtPlot or rtOTF, and they don't use this routine.
%
% Examples:
%   [otf,fSupport] = oiCalculateOTF(oi);
%   [otf,fsmm] = oiCalculateOTF(oi,'wave','mm');
%
% See also:  applyOTF, oiCalculateOTF, oiCompute, dlMTF, dlCore
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('wave'), wave = oiGet(oi, 'wave'); end
if notDefined('unit'), unit = 'cyclesPerDegree'; end

optics = oiGet(oi,'optics');
opticsModel = opticsGet(optics,'model');

% Retrieve the frequency support in the proper units.
fSupport = oiGet(oi, 'frequencysupport',unit);

switch lower(opticsModel)
    case {'dlmtf','diffractionlimited'}
        otf = dlMTF(oi,fSupport,wave,unit);
    
    case {'custom','shiftinvariant'}
        % Calculate the OTF at each wavelength from the custom data.
        % Also return these OTFs at the specified frequency
        % support, which depends on the optical image.
        %
        % It is important that the units specified for this calculation and
        % the units specified for the custom OTF be the same
        otf = customOTF(oi, fSupport, wave, unit);

    case {'skip','skipotf'}
        % Doesn't really happen.
        warning('No OTF method selected.');
    
    otherwise
        error('Unknown optics model: %s',opticsModel);
end

end