function OTF2D = dlCore(rho, inCutFreq)
% Compute 2D diffraction limited optical transfer function
%
% Syntax:
%   OTF2D = dlCore(rho, inCutFreq)
%
% Description:
%    This routine calculates the 2D optical transfer function (OTF) of a
%    diffraction limited optical system. The OTF is returned with DC in the
%    (1, 1) (upper left) position. Apply fftshift to put DC in the center.
%
%    The 2D OTF depends only a few parameters. These are the radial
%    distance (rho) of each spatial frequency from the origin, and the
%    input cutoff frequency (inCutFreq) at which wavelength. The inCutFreq,
%    in turn, only depends on the f# of the diffraction limited system.
%
%    The variable rho is a matrix defining the frequency of each entry
%    (cycles per meter). The vector inCutFreq is the incoherent cutoff
%    frequency as a function of wavelength.
%
%    A little more formally, an aberration-free diffraction-limited system
%    with a circular exit pupil can be described as:
%
%       di = distance between aperture and detector (meters)
%       A = aperture diameter  (meters)
%       lambda = wavelength of incident light (meters)
%       rho = frequency (cycles per meter)
%       rho0 = (A / 2 * lambda * di) (cycles/meter)
%
%    The formula for the OTF at frequency rho and wavelength lambda is
%    often quoted as
%
%     H(rho, lambda)
%        = (2 / pi) * (acos(rho / (2 * rho0)) - ...
%           (rho / 2 * rho0) * sqrt(1 - (rho / (2 * rho0) ^ 2)))
%        = (2 / pi) * (acos(rho / inCutF)   - ...
%           (rho / inCutF) * sqrt(1 - (rho / inCutF)) ^ 2))
%        or 0 if rho >= 2 * rho0, which is rho / inCutFreq >= 1
%
%    Simplify using the spatial cutoff frequency (2 * rho0)
%
%      inCutFreq = (A / (di * wavelength)) = 2 * rho0
%
%    Define the normalized frequency, rho/inCutF.
%
%       normF = rho / (A / (d i *wavelength))
%
%    In that case, the OTF formula becomes
%
%       H(normF, lambda) = (2 / pi) * ...
%           (acos(normF) - normF * sqrt((1 - normF) ^ 2))
%
%    In this form, we convert normalized frequency OTF back to real units
%
% Inputs:
%    rho       - Numeric. The radial distance of each frequency from the
%                origin, in either vector or scalar form.
%    inCutFreq - Numeric. Input cutoff frquency wavelength
%
% Outputs:
%    OTF2D     - Numeric. A Two-Dimensional OTF, of a size dependent on the
%                input's scale. A scalar rho and inCutFreq result in a
%                scalar OTF2D, whereas a vector for rho results in a vector
%                for the OTF2D. If both variables are not scalar, a matrix
%                is returned.
%
% Optional key/value pairs:
%    None.
%
% Sources:
%    * http://ao.osa.org/ViewMedia.cfm?id=38173&seq=0 Muralidhara Subbarao,
%      APPLIED OPTICS / Vol. 29, No. 4 / 1 February 1990 Optical transfer
%      function of a diffraction-limited system for polychromatic
%      illumination.
%    * http://www.microscopyu.com/articles/optics/mtfintro.html
%
% See Also:
%   dlMTF
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/08/18  jnm  Formatting
%    06/28/19  JNM  Documentation update

% Number of spatial frequency and wavelength samples
[r, c] = size(rho);
nWaves = length(inCutFreq);

% Pre-allocate the OTF2D.
OTF2D = zeros(r, c, nWaves);

% The OTF formula used here places DC in the center of the matrix. We want
% DC to be in the (1, 1) position, so we apply an ifftshift at the last
% step This part can be boosted by gpu computing
for ii = 1 : nWaves
    nFreq = rho / inCutFreq(ii);
    nFreq(nFreq > 1) = 1;
    otf = (2 / pi) * (acos(nFreq) - nFreq .* sqrt(1 - nFreq .^ 2));
    OTF2D(:, :, ii) = ifftshift(otf);
end

% If there is only one wavelength, return the data as a single image.
if nWaves == 1, OTF2D = squeeze(OTF2D); end

end