function [OTF2D, fSupport] = customOTF(oi, fSupport, wavelength, units)
% Interpolate optics OTF for shift-invariant calculation in optical image
%
% Syntax:
%   [OTF2D, fSupport] = customOTF(oi, [fSupport], [wavelength], [units])
%
% Description:
%    In the shift-invariant optics model, custom data are stored in the
%    optics.OTF slot. This routine reads the otf data and interpolates them
%    to the fSupport and wavelength of the optical image structure or
%    optics structure.
%
%    The returned OTF is normalized so that all energy is transmitted
%    (i.e., the DC value is 1). This is done by normalizing the peak value
%    to one. If we ever have a case when the peak is other than the DC, we
%    have a problem with energy conservation - where did the photons go?
%
%    The units for the frequency support are cycles/millimeters.
%    Perhaps we should add a 'units' input argument here.
%
% Inputs:
%    oi         - Struct. An optical image structure
%    fSupport   - (Optional) Frequency support. Default retrieves from oi
%    wavelength - (Optional) Vector. Wavelengths. Default retrieves from oi
%    units      - (Optional) String. OTF units. Default millimeter ('mm')
%
% Outputs:
%    OTF2D      - Matrix. Two-dimensional OTF for each wavelength.
%    fSupport   - Matrix. Frequency support.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   oiCalculateOTF, dlMTF, dlCore
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/29/17  dhb  Comment cleaning. Replace outer fftshift around
%                   interpolation with ifftshif, in the branch that loops
%                   over wavelength. ifftshift is what we want here, and
%                   indeed that was correctly used in the branch that
%                   didn't loop over wavelength.
%    03/08/18  jnm  Formatting
%    06/28/19  JNM  Documentation update

% Handle optional args
if notDefined('oi'), error('Optical image required.'); end
if notDefined('wavelength'), wavelength = oiGet(oi, 'wavelength'); end

% In the custom case, we think the units should always be millimeters.
if notDefined('units'), units = 'mm'; end
if notDefined('fSupport'), fSupport = oiGet(oi, 'fSupport', units); end

% Get the frequency support
fx = fSupport(:, :, 1);
fy = fSupport(:, :, 2);
nX = size(fx, 2);
nY = size(fy, 1);

% Get number of wavelengths
nWave = length(wavelength);

% Grab the optics structure from the oi
optics = oiGet(oi, 'optics');

% Get frequency representaton of OTF in the optics structure
%
% The units should be specified and optional; we are using mm for now.
otfSupport = opticsGet(optics, 'otfSupport');
[X, Y] = meshgrid(otfSupport{1}, otfSupport{2});

% Set interpolation method
if ieSessionGet('gpu compute')
    method = 'linear';
    fx = gpuArray(fx);
    fy = gpuArray(fy);
else
    method = '*linear';
end

% Find OTF at each wavelength. We may be interpolating from the custom data
if isscalar(wavelength)
    % Not entirely clear we need to interpolate in this case, as the
    OTF2D = opticsGet(optics, 'otfData', wavelength);

    % Do interpolation, use fftshift/ifftshift to take care of DC
    % positions. The reason we have to do this is that the frequency
    % support is represented with DC in the center, while the OTF is stored
    % with DC at (1, 1). That is, the frequency support in the optics
    % object matches the OTF only after application of fftshift.
    %
    % The 0 as the last argument of the interp2 call is the value used for
    % extrapolation. Putting this here prevents NaNs should the requested
    % support exceed the input range
    OTF2D = ifftshift(interp2(X, Y, fftshift(OTF2D), fx, fy, method, 0));
else
    % Do it wavelength by wavelength
    if ieSessionGet('gpu compute')
        OTF2D = zeros(nY, nX, nWave, 'gpuArray');
    else
        OTF2D = zeros(nY, nX, nWave);
    end
    for ii = 1:length(wavelength)
        tmp = opticsGet(optics, 'otfData', wavelength(ii));

        % See comments above for how this interpolation works.
        OTF2D(:, :, ii) = ...
            ifftshift(interp2(X, Y, fftshift(tmp), fx, fy, method, 0));
    end
end

end