function optics = opticsMouse(pupilRadius)
% Optics for mouse eye
%
% Syntax:
%   optics = opticsMouse(pupilRadius)
%
% Description:
%    Not yet tested or working. Just a placeholder for now.
%
%    Optics calculated for a mouse's eye. This is similar to human optics,
%    with different parameter values.
%
%    There are extensive notes contained in this function. To access the
%    notes type 'edit opticsMouse.m' into the Command Window.
%
% Inputs:
%    pupilRadius - Numeric. The pupil radius in meters.
%
% Outputs:
%    optics      - Structure. The optics structure for the mouse eye.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16       ISETBIO Team, 2016
%    03/12/18  jnm  Formatting

%%
if notDefined('pupilRadius'), pupilRadius = 0.00059; end
fLength = 0.001756;  % 1.756 mm
%% NOTES ABOUT FOCAL LENGTH
% Since the mouse optical system is not symmetric, there are two focal
% lengths: front (before the cornea + lens) and back (after the cornea +
% lens, roughly the lens-retina distance).
%
% Since light travels in air before the optical system, and in vitreous
% humor after, the value for the back focal length should be multiplied by
% the refraction index to be used. We have:
%
%       OpticalPower = 1/focalLength*n_medium
%                    = 1/frontFocalLength*n_air = 1/0.001756*1
%                    = 1/backFocalLength*n_vitreous = 1/0.002347*1.3341
%
% Note: this value of n_vitreous = 1.3341 comes from geometrical
% calculations on the mouse eye. For comparison, the human vitreous humor
% has an index of 1.336.
%
% Therefore, to avoid having to multiply all the geometry by the refraction
% index, we use the front focal length, for which the refraction index is
% 1.0003, approx to 1.
%
% The optical power is wavelength-dependent, which means the focal length
% is too: longer focal length for short wavelengths, shorter focal length
% for long wavelengths. For fLength, we need the focal length at the
% infocus wavelength.
% (DO WE?? TODO: what is the dioptric power of the unaccomodated eye?)
%
% (i.e. that focuses on the retina surface, or as close to it as possible.)
% At 544nm, we have 0.4 diopters of defocus, which is our best value. (see
% the plot in mouseCore.m for defocus values. By interpolation, the zero
% defocus would be at about 542nm.)
%
% So ***we are using the front focal length at 544nm for fLength value***.
%
% (Source: A SCHEMATIC EYE FOR THE MOUSE, AND COMPARISONS WITH THE RAT,
% Remtulla and Hallett, 1984)

optics.type = 'optics';
optics.name = 'mouse';
optics = opticsSet(optics, 'model', 'shiftInvariant');

% Ratio of focal length to diameter.
optics = opticsSet(optics, 'fnumber', fLength / (2 * pupilRadius));
optics = opticsSet(optics, 'focalLength', fLength);

% The OTF computation is the same as the human, except the values for
% dioptricPower and pupilRadius.
optics = opticsSet(optics, 'otfMethod', 'mouseOTF');
% Compute the OTF and store it. We use a default pupil radius, dioptric
% power, and so forth.
% diopter power: About 570 diopters for the mouse, 60 for the human
dioptricPower = 1 / fLength;

% Try to get current wavelengths: in the current oi, or the current scene.
oi = vcGetObject('oi');
if isempty(oi)
     getScene = 1;
else
    spect = oiGet(oi, 'wave');
    if isempty(spect)
         getScene = 1;
    else
        wave = spect;
        getScene = 0;
    end
end
if getScene
    scene = vcGetObject('scene');
    if isempty(scene)
        wave = (325:5:635)';
    else
        spect = sceneGet(scene, 'wave');
        if isempty(spect)
            wave = (325:5:635)';
        else
            wave = spect;
        end
    end
end
optics = opticsSet (optics, 'wave', wave);

%% The human optics are an SI case, and we store the OTF at this point.
[OTF2D, frequencySupport] = mouseOTF(pupilRadius, dioptricPower, [], wave);
optics = opticsSet(optics, 'otfData', OTF2D);

% [Note: XXX - This OTF is a diffraction-limited model, with mouse values
% for defocus. I'm not sure what happens for the chromatic aberration? The
% Williams factor (used in humanOTF) is skipped for the mouseOTF.]
% Another possibility would be to use a diffraction-limited model, with the
% mouse values for pupil-size and focal length, and add chromatic
% aberration from data (??), and build it all with Zmax. Peter Catrysse
% said he could probably help with that...

%% Frequency is in cyc/mm
metersPerDegree = fLength * tan(1 / 180 * pi);
millimetersPerDegree = metersPerDegree * 1000;
% The frequencies are in cycles/degree. We convert to cycles/mm here.
frequencySupport = frequencySupport * (1/millimetersPerDegree);
fx = frequencySupport(1, :, 1);
fy = frequencySupport(:, 1, 2);
optics = opticsSet(optics, 'otffx', fx(:)');
optics = opticsSet(optics, 'otffy', fy(:)');

optics = opticsSet(optics, 'otfWave', wave);

% Transmittance of the lens: get it from data file, interpolate to fit
% wavelengths if needed
fname = '~/psych221/mouseTransmittance.mat';
extrapVal = 0; % transmittance is 0 outside of the known mouse wavelength
res = ieReadSpectra(fname, wave, extrapVal);
optics = opticsSet(optics, 'transmittance', res); % real transmittance
% no transmittance: same on all wavelengths
%optics = opticsSet(optics, 'transmittance', ones(length(wave), 1));

% figure(1);
% mesh(frequencySupport(:, :, 1), frequencySupport(:, :, 2), ...
%    OTF2D(:, :, 20));
% mesh(abs(otf2psf(OTF2D(:, :, 15))))
%

end
