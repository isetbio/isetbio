% Illustrate photoPigment object
%
% Description:
%    The photoPigment object represents the spectral responsivity of
%    the cone photopigments.  These values are needed to calculate the
%    capture of light by cones, along with the oi spectral irradiance
%    and the other inert pigments (lens, macular). 
%
%    This tutorial illustrates how it works.
%
% See also
%    

%% Initialize
ieInit; clear; close all;

%% Set wavelength support
wls = (400:5:700)';

%% Generate the ring rays stimulus
radialFrequency = 8;
pixelSize = 256;
scene = sceneCreate('rings rays',radialFrequency,pixelSize,wls);
scene = sceneSet(scene, 'fov', 0.5);

%% Compute the retinal image
%
% First create the optical image object.
pupilDiamterMm = 3;
%oi = oiCreate('wvf human', pupilDiamterMm , [], wls);
oi = oiCreate('wvf human', pupilDiamterMm);

% Compute the optical image
oi = oiCompute(scene, oi);

% The lens density in the oi object affect the cone fundamentals.  This
% code shows how to adjust lens density.

% The values in the unit density variable are multiplied by the scalar in
% the density variable to get the overall spectral density, and from thence the
% transmittance.  See "help Lens".
lens0 = oiGet(oi,'lens');
lensUnitDensity1 = lens0.unitDensity;
lensPeakDensity1 = lens0.density;
lens1 = Lens('wave',wls,'unitDensity',lensUnitDensity1,'density',lensPeakDensity1);
oi = oiSet(oi,'lens',lens1);



%% Generate the mosaic
cm = cMosaic(...
    'sizeDegs', [1.0 1.5], ...     % SIZE: 1.0 degs (x) 0.5 degs (y)
    'eccentricityDegs', [0 0], ...  % ECC: (0,0)
    'eccVaryingConeBlur', true ...
    );

%% Create photoPigment object
% Plot spacing at 5 nm looks a little nicer than the default spacing.
pp = photoPigment('wave', 400:5:700);

%% Wavelength
% The internal wavelength representation of the parameters is from 390:830
% at 1 nm. The user can set the wavelength representation used for
% interface with other routines. By default this is 400:10:700.

%% Absorbance
% The cone absorbance function used by default is obtained via data routine
% coneAbsorbanceReadData.
%
% Absorbance is sometimes called optical density.
%
% The peak absorbance is 1 by convention in how the values are tabulated.
ieNewGraphWin;
plot(pp.wave, pp.absorbance)
grid on;
xlabel('Wavelength (nm)');
ylabel('Relative sensitivity');

%% Geometry
%
% See coneDensityReadData for explanation of retinal coordinate system.
eccentricity = 0.0;
angle = 0;

% Spacing (um), aperture (um), density (cones/mm2)
[s, a, d] = coneSizeReadData('eccentricity', eccentricity, 'angle', angle);

%%