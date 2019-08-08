function [OTF2D, fSupport, inCutoffFreq] = dlMTF(oi, fSupport, wave, units)
% Gateway routine that assembles params to compute diffraction limited OTF
%
% Syntax:
%   [OTF2D, fSupport, inCutoffFreq] = ...
%       dlMTF(oi, [fSupport], [wave], [units])
%
% Description:
%    Compute the diffraction limited 2D OTF (OTF2D) at each wavelength for
%    the optics within the optical image. The OTF is returned with DC in
%    the (1, 1) position. This is the upper left. Apply fftshift to put DC
%    in the center.
%
%    The diffraction limited OTF only depends on the optics. But for some
%    units and conditions we need to know properties of the optical image
%    to perform the calculation. If you know these parameters, we make it
%    possible to call this routine using dlMTF(optics, ...). This is
%    possible because we check at the beginning of the routine to see
%    whether the first argument is of type optical image or of type optics.
%
%    The frequency support(fSupport) and the incoherent cutoff frequency as
%    a function of wavelength (in nm) can also be calculated and returned.
%    The units for the frequency support, cycles/{meters, millimeters,
%    microns}, can be specified (units).
%
%    The formulae are described in dlCore.m
%
%    This function contains examples of usage inline. To access these, type
%    'edit dlMTF.m' into the Command Window.
%
% Inputs:
%    oi           - Struct. An optical image structure.
%    fSupport     - (Optional) Matrix. The fSupport function. Default is to
%                   retrieve the oi's fSupport with the provided units.
%    wave         - (Optional) Vector. A vector containing the wavelengths.
%                   Default is retrieve from the provided oi.
%    units        - (Optional) String. A string describing the units.
%                   Default is cyclesPerDegree.
% Outputs:
%    OTF2D        - Matrix. A 3D Matrix containing the 2D diffraction
%                   limited OTF at each wavelength for the provided oi.
%    fSupport     - Matrix. The fSupport function for the OTF.
%    inCutoffFreq - Vector. A cutoff frequency vector for each wavelength.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   oitCalculateOTF, dlCore.m
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/08/18  jnm  Formatting
%    04/07/18  dhb  Deep six broken examples.
%    06/28/19  JNM  Documentation update

% Examples:
%{
    oi = oiCreate('diffraction limited');
    wavelength = 700;
    unit = 'mm';
    fSupport = oiGet(oi, 'fSupport', unit);
    OTF2D = dlMTF(oi, fSupport, wavelength, unit);
    vcNewGraphWin;
    mesh(fSupport(1, :, 1), fSupport(:, 1, 2), fftshift(OTF2D));
    colorbar;
    xlabel('cyc/mm');
    ylabel('cyc/mm');
%}
%{
    scene = sceneCreate;
    oi = oiCreate('diffraction limited');
    oi = oiCompute(oi, scene);
    OTF2D = dlMTF(oi);
    size(OTF2D)
    vcNewGraphWin;
    mesh(fftshift(OTF2D(:, :, 1)));
    vcNewGraphWin;
    mesh(fftshift(OTF2D(:, :, end)));
%}

%%
if notDefined('oi'), error('Optics or optical image required.'); end

% The user can send in the OI or OPTICS. We only need OPTICS.
if strcmpi(opticsGet(oi, 'type'), 'opticalimage')
    optics = oiGet(oi, 'optics');
else
    if nargin < 4
        error('All dlMTF arguments must be sent if optics is first');
    end
    optics = oi;
    clear oi
end

if notDefined('wave'), wave = oiGet(oi, 'wavelength'); end
if notDefined('units'), units = 'cyclesPerDegree'; end
if notDefined('fSupport'), fSupport = oiGet(oi, 'fSupport', units); end

apertureDiameter = opticsGet(optics, 'aperture diameter');
fpDistance = opticsGet(optics, 'focal Plane Distance');

fx = fSupport(:, :, 1);
fy = fSupport(:, :, 2);

%  Distance of each frequency pair from the origin (rho, theta)
%  (dc = [0, 0]). The frequency support is in cycles/deg.
rho = sqrt(fx .^ 2 + fy .^ 2);

% Wavelength is stored in nanometers. This converts it to meters, the same
% units as the apertureDimaeter.
wave = wave *  1e-9;

% This formula assumes that a lens that is free of spherical aberration and
% coma. When the source is far away, fpDistance is the focal length and
% the ratio (apertureDiameter / fpDistance) is the f#.
%
%  see discussion in dlCore.m
%
% http://spie.org/x34304.xml - Cutoff frequency
% http://spie.org/x34468.xml - Airy disk
%
inCutoffFreq = (apertureDiameter / fpDistance) ./ wave;

switch lower(units)
    case {'cyclesperdegree', 'cycperdeg'}
        % cycle/meter * meter/deg -> cycles/deg
        % Used for human calculations?
        inCutoffFreq = inCutoffFreq * ...
                    oiGet(oi, 'distancePerDegree', 'meters');

    case {'meters', 'm', 'millimeters', 'mm', 'microns', 'um'}
        inCutoffFreq = inCutoffFreq/ieUnitScaleFactor(units);

    otherwise
        error('Unknown units.');
end

% Now, both rho and the cutoff frequency are in cycles/degree.
OTF2D = dlCore(rho, inCutoffFreq);

end