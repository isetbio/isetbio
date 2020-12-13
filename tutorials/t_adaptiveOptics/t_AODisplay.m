%% Show how to create a display that models some AO experiments
%
% Description:
%

% History:
%   12/03/20  dhb  Added comments.

%% Initialize
close all; clear; ieInit;

%% Set up wavelength support for stimulus
wls = (700:900)';
nWls = length(wls);

% Stimulus spaital parameters
backgroundSizeDegs = 0.85;                  % Linear side of square background field.  This can be the imaging field.
testDiameterDegs = 0.1;                     % Diameter of circular test spot. 
nPixels = 512;                              % Number of h and v pixels in background

% Stimulus wavelengths
backgroundWlNm = 840;                       % Wavelength of monochromatic background field. This can be the imaging light
testWlNm = 830;                             % Wavelenght of monochromatic test spot.

%% Stimulus power.
%
% Specified here as the power entering the eye (going through
% the pupil).  If the beam just fills the pupil, then this is just the power
% measured at the cornea.
%
% If you have instead corneal irradiance in UW/mm^2 and this fills or overfills
% the pupil, then multiply by pupil area in mm^2 to get the power entering the eye. 
%
% This code assumes that the power for both background is spread out over
% the background area. This is because typically we measure power for the
% whole raster, and then in an experiment switch it off for part of the
% time to make the spatial pattern we want.
%
% The test should be specified as the amount of power added to
% the background (its incremental power).  This is what you'll get if
% background and test are in separate channels of the system and you
% measure them separately.
pupilDiameterMm = 7;                            % Pupil size. 
backgroundCornealPowerUW = 10;              % Power passing through the pupil, measured at the cornea, for background light.
testCornealPowerUW = 10000;                 % Power passing through the pupil, measured at the cornea, for the test spot

% Get equivalent spectral radiance of background and test increment as a 
% function of the wavelength support.
%
% The routine here finds the radiance on an external conventional display that
% produces the same retinal illuminance as the corneal power specified
% above.  This is purely geometric calculation; attenuation of light by
% occular media is not taken into account at this stage.
%
% If you already have the stimulus spd as radiance in units of
% Watts/[sr-m2-nm], you can skip this call and simply set those functions
bgRadiance = AOMonochromaticCornealPowerToRadiance(wls,backgroundWlNm,backgroundCornealPowerUW,pupilDiameterMm,backgroundSizeDegs^2);
testRadiance = AOMonochromaticCornealPowerToRadiance(wls,testWlNm,testCornealPowerUW,pupilDiameterMm,backgroundSizeDegs^2);

% Create an empty scene to use for the spot.  Put it far enough away so it
% is in focus for an emmotropic eye accommodated to infinity.
theScene = sceneCreate('empty');
theScene = sceneSet(theScene,'wavelength',wls);
theScene =  sceneSet(theScene,'distance',2);
theScene = sceneSet(theScene, 'h fov', backgroundSizeDegs);

% Make an image with the background + spot spectral radiance at all locations
%
% First use extant routine drawSpot to draw the spot. 
spatialParams.row = nPixels;
spatialParams.col = nPixels;
spatialParams.type = 'Spatial';
spatialParams.spatialType = 'spot';
spatialParams.fieldOfViewDegs = backgroundSizeDegs;
spatialParams.backgroundSizeDegs = backgroundSizeDegs;
spatialParams.spotSizeDegs = testDiameterDegs;
spotPattern = drawSpot(spatialParams);

% Then fill in appropriate radiance at each pixel
radianceEnergySpot = zeros(spatialParams.row,spatialParams.col,nWls);
for i = 1:spatialParams.row
    for j = 1:spatialParams.col
        % Background pixels are 1, spot pixels are 2
        if (spotPattern(i,j) == 1)
            radianceEnergySpot(i,j,:) = bgRadiance;
        elseif (spotPattern(i,j) == 2)
            radianceEnergySpot(i,j,:) = bgRadiance + testRadiance;
        end
    end
end

% Convert radiance to quantal units
radiancePhotonsSpot = Energy2Quanta(wls,radianceEnergySpot);

% Put in the image in quantal units (photons)
theScene = sceneSet(theScene,'photons',radiancePhotonsSpot);

% Visualize scene
vcAddAndSelectObject(theScene);
sceneWindow;

%% Build a simple retinal image
defocusAmount = 0.1;
pupilDiameterMm = 7;
accommodatedWavelength = 700;
zCoeffs = zeros(66,1);
defeatLCA = false;

% Set up wavefront optics object
%
% Compute pupil function using 'no lca' key/value pair to turn off LCA.
% You can turn it back on to compare the effect.
%
% Deal with best focus by specifying that the wavefront parameters
% were measured at the wavelength we want to say is in focus. This
% is a little bit of a hack but seems OK for the diffraction limited case
% we're using here.
wvfP = wvfCreate('calc wavelengths', wls, 'zcoeffs', zCoeffs, ...
    'name', sprintf('human-%d', pupilDiameterMm));
wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMm);
wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWavelength);
wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMm);
wvfP = wvfSet(wvfP,'zcoeffs', defocusAmount, 'defocus');

% Compute pupil function and PSF
%
% Whether LCA should be included depends on your apparatus and
% is controlled by the boolean defeatLCA in the computation of
% the pupil function.
wvfP = wvfComputePupilFunction(wvfP,false,'no lca',defeatLCA);
wvfP = wvfComputePSF(wvfP);

% Generate optical image object from the wavefront object
theOI = wvf2oi(wvfP);

% Set some properties of the object, mainly to show how we
% extract, set, and put this back.
opticsP = oiGet(theOI, 'optics');
opticsP = opticsSet(opticsP, 'model', 'shift invariant');
opticsP = opticsSet(opticsP, 'name', 'human-wvf-nolca');
theOI = oiSet(theOI,'optics',opticsP);

% wave = [450:100:950]';
% % theOI = oiCreate('wvf human');
% wvfP = wvfCreate('calc wavelengths', wave, 'zcoeffs', zCoeffs, 'measured pupil', pupilDiamMm,'calc pupil', pupilDiamMm);
% wvfP = wvfSet(wvfP,'zcoeffs', defocusAmount, 'defocus');
% wvfP = wvfComputePSF(wvfP);
% oi_AOSLO = wvf2oi(wvfP);

% Compute the retinal image and visualize
theOI = oiCompute(theOI, theScene);
vcAddAndSelectObject(theOI);
oiWindow;

%% Create the coneMosaic object
% Generate default human optics 
%defaultHumanOI = oiCreate('wvf human');

% Generate a human foveal cone mosaic
resamplingFactor = 7;
integrationTimeSecs = 50/1000;
cMosaic = coneMosaicHex(resamplingFactor , ...
    'wave',wls, ...
    'eccBasedConeDensity', false, ...
    'eccBasedConeQuantalEfficiency', false, ...
    'maxGridAdjustmentIterations',200, ...
    'integrationTime', integrationTimeSecs, ...
    'fovDegs', [1.2*backgroundSizeDegs 1.2*backgroundSizeDegs]);

% Eye movement path.  Here, just one time point and no motion.
nTrialsNum = 1;
emPath = zeros(nTrialsNum, 1, 2);

% Compute isomerizations for each eye position.
cMosaic.compute(theOI);
cMosaic.window;