%% Show how to create a display that models some AO experiments
%
% Description:
%    This tutorial shows how to set up an ISETBio scene that represents a
%    display used in an adaptive optics apparatus. It has various options that 
%    illustrate how different displays might be simulated.  It computes
%    cone isomerizations for a spot against a background.
%
%    In particular note options to vary stimulus spot wavelength, to
%    specify stimulus intensity in terms of conrneal power or radiance of a
%    display, to turn modeling of longitudinal chromatic aberration on or
%    off, and to provide custom specification of paramters such as lens
%    density, macular pigment density, and cone photopigment absorbance.
%    The latter may requires some substantive thought for accurate modeling
%    of what happens out at the near IR wavelengths used in AO experiments.
%    Although standard tables go out to about 830 nm, it is not immediately
%    clear how accurate these tables are at such long wavelengths.
%       
%    This tutorial has not worried much about how realistic these
%    specified powers are, nor whether they are within light safety
%    limits. Nor has it been highly tested.  But the basic operation seems
%    in place
%
% See also: t_noLcaOptics, t_radianceToCornealIrradiance, cMosaicHex, Lens,
%           Macular, receptorPigment

% History:
%   12/03/20  dhb  Added comments.
%   12/27/20  dhb  Lots of tune-up and comments.
%   04/11/21  dhb  Correctly set pupil size in optical image structure.
%                  This was not happening, because it isn't done in
%                  wvf2oi().

%% Initialize
close all; clear; ieInit;

%% Stimulus spectral properties
%
% Narrowband stimulus modeled as Gaussians centered on specified
% wavelength, with specified FWHM
bgWl = 840;                               % Wavelength of narrowband background field. This can be the imaging light
bgFWHM = 21;                              % Full-width at half max of background (in nm).

% Choices of test wavelength currently are:
%   550
%   840
% It would be easy to add more, the limit is just so I could hand set
% reasonablish intensities and wavelength support for each choice.
testWl = 550;                             % Wavelength of monochromatic test spot.
testFWHM = 26;                            % Full-width at half max of test (in nm).

%% Set up wavelength support for stimulus
% 
% Wavelength spacing set to 5 to speed up calculations
% a bit.  Spot factor handles test wavelength specific
% intensity of spot relative to background.
deltaWl = 5;
switch (testWl)
    case (550)
        wls = (500:deltaWl:900)';
        spotFactor = 1/50;
    case (840)
        wls = (700:deltaWl:900)';
        spotFactor = 5;
    otherwise
        error('Need to add test wavelength to switch statements in tutorial')
end
nWls = length(wls);

%% Stimulus spatial parameters
backgroundSizeDegs = 0.85;                  % Linear side of square background field.  This can be the imaging field.
testDiameterDegs = 0.1;                     % Diameter of circular test spot.
nPixels = 512;                              % Number of h and v pixels in background field (aka full field)

% Make relative spectral power distributions. Each is approximated by a
% Gaussian with specified center wavlength and FWHM, and with total power
% given by the corneal power specified above.  The call to trapz takes
% wavelength spacing into account when normalizing the power.
bgRelSpd = normpdf(wls,bgWl,FWHMToStd(bgFWHM));
bgUnitSpd = bgRelSpd/trapz(wls,bgRelSpd);
testRelSpd = normpdf(wls,testWl,FWHMToStd(testFWHM));
testUnitSpd = testRelSpd/trapz(wls,testRelSpd);

%% Pupil size
pupilDiameterMm = 7;                        % Pupil diameter.
pupilAreaMm = pi*(pupilDiameterMm/2)^2;     % For convenience below, compute pupil area.

%% Stimulus power
%
% How are we specifing power?
% Options are:
%   'cornealPower' - Specify power passing through the pupil in UW.
%   'radiance'     - Specify the radiance of a display
powerSpecification = 'radiance';

% Always specify background corneal power, just to keep other choices
% locked to it by appropriate precomputed constants.  See comment below for
% convention on this specification.
bgCornealPowerUW = 10;

% Other aspects get handled through this switch depending on what we are
% demonstrating.
switch (powerSpecification)
    case 'cornealPower'
        % Specified here as the power entering the eye (going through
        % the pupil).  If the beam just fills the pupil, then this is just the power
        % measured at the cornea.
        %
        % If you have instead corneal irradiance in UW/mm^2 and this fills or overfills
        % the pupil, then multiply by pupil area in mm^2 to get the power entering the eye.
        %
        % This code assumes that the power for both background and test is spread
        % out over the background area. This is because typically for an AOSLO we
        % measure power for the whole raster, and then in an experiment switch it
        % off for part of the time to make the spatial pattern we want.
        %
        % The test should be specified as the amount of power added to
        % the background (its incremental power).  This is what you'll get if
        % background and test are in separate channels of the system and you
        % measure them separately.
        testCornealPowerUW = spotFactor*bgCornealPowerUW;   % Power passing through the pupil, measured at the cornea, for the test spot.
           
        % Make spds that give desired full power.
        backgroundSpdCornealPowerUW = bgCornealPowerUW*bgUnitSpd;
        testSpdCornealPowerUW = testCornealPowerUW*testUnitSpd;
        
        % Get equivalent spectral radiance of background and test increment as a
        % function of the wavelength support.
        %
        % The routine here finds the radiance on an external conventional display that
        % produces the same retinal illuminance as the corneal power specified
        % above.  This is purely geometric calculation; attenuation of light by
        % occular media is not taken into account at this stage.  Note that this
        % conversion routine expects power per wavelength band, not power per nm,
        % as its input, but returns power per nm as its output.  A little
        % confusing, but that is why the spd being passed in is multiplied by the
        % wavelength spacing.
        bgSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,backgroundSpdCornealPowerUW*deltaWl,pupilDiameterMm,backgroundSizeDegs^2);
        testSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,testSpdCornealPowerUW*deltaWl,pupilDiameterMm,backgroundSizeDegs^2);
        
        % Make sure our computed radiance yields the desired corneal
        % irradiance when we go in the other direction.
        bgSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(bgSpdRadiance,backgroundSizeDegs^2)*(1e6)*((1e-3)^2);
        bgCornealPowerUWCheck = trapz(wls,bgSpdCornealIrradianceUWMm2Check)*pupilAreaMm;
        if (abs(bgCornealPowerUWCheck-bgCornealPowerUW)/bgCornealPowerUW > 1e-4)
            error('Do not get right cornal power back from computed radiance');
        end
        testSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(testSpdRadiance,backgroundSizeDegs^2)*(1e6)*((1e-3)^2);
        testCornealPowerUWCheck = trapz(wls,testSpdCornealIrradianceUWMm2Check)*pupilAreaMm;
        if (abs(testCornealPowerUWCheck-testCornealPowerUW)/testCornealPowerUW > 1e-4)
            error('Do not get right cornal power back from computed radiance');
        end
        
    case 'radiance'
        % Provide the total radiance in Watts/[sr-m2] of the stimulus being
        % viewed.  The numbers here correspond to the cornal powers
        % specified above, for the default geometry.
        bgRadiance = bgCornealPowerUW*(1.1807e+03/10);
        testRadiance = spotFactor*bgRadiance;
        
        % Convert to full spectral radiance 
        bgSpdRadiance = bgRadiance*bgUnitSpd;
        testSpdRadiance = testRadiance*testUnitSpd;
    otherwise
        error('Unknown power specification given');
end

%% Create scene
%
% Create an empty scene to use for the spot.  Put it far enough away so it
% is basically in focus for an emmetropic eye accommodated to infinity.
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
            radianceEnergySpot(i,j,:) = bgSpdRadiance;
        elseif (spotPattern(i,j) == 2)
            radianceEnergySpot(i,j,:) = bgSpdRadiance + testSpdRadiance;
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
%
% Assume observer is accommodated to the test wavelength.
% You have the option here of turning off LCA, if that is the
% way your AO system is set up.  It is left on by default.
defocusAmount = 0.1;
pupilDiameterMm = 7;
accommodatedWl = testWl;
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
wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWl);
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

% Set the fNumber to correspond to the pupil size
focalLengthMM = oiGet(theOI,'focal length')*1000;
theOI = oiSet(theOI, 'optics fnumber', focalLengthMM/pupilDiameterMm);

% You could adjust lens density if you don't like the default values for
% the very long wavelengths used here.  But, the default density at long
% wavelengths is zero, which is probably what you want. Change the values
% in variables lensUnitDensity1 and lensPeakDensity1 to do this to do this,
% as illustrated by the conditional.
%
% The values in the unit density variable are multiplied by the scalar in
% the density variable to get the overall spectral density, and from thence the
% transmittance.  See "help Lens".
%
% Note that this code is not highly tested, but it runs without crashing.
lens0 = oiGet(theOI,'lens');
if (any(lens0.wave ~= wls))
    error('Wavelengths in Lens object are not as expected');
end
lensUnitDensity1 = lens0.unitDensity;
lensPeakDensity1 = lens0.density;

% Set to true to see that making lens opaque at the test wavelength
% attenuates the test spot.  The effect will be clearest in the
% isomerizations image, because the optical image display is normalized and
% will show the remaining spot from the tails of the Gaussian spd.
changeLens = false;
if (changeLens)
    wlsChange = find(testSpdRadiance > 0.001*max(testSpdRadiance));
    lensUnitDensity1 = zeros(size(lensUnitDensity1));
    lensUnitDensity1(wlsChange) = 1;
    lensPeakDensity1 = 100;
end
lens1 = Lens('wave',wls,'unitDensity',lensUnitDensity1,'density',lensPeakDensity1);
theOI = oiSet(theOI,'lens',lens1);

% Set some properties of the object, mainly to show how we
% extract, set, and put this back.
opticsP = oiGet(theOI, 'optics');
opticsP = opticsSet(opticsP, 'model', 'shift invariant');
opticsP = opticsSet(opticsP, 'name', 'human-wvf-nolca');
theOI = oiSet(theOI,'optics',opticsP);

% Compute the retinal image and visualize.  Note that to convolve with the
% PSF, we have to pad the outside of the image with something, and we
% haven't fussed here with the defaults in how that's done.  That isn't a
% problem as long as you recognize that you shouldn't interpret the edges
% of the computed retinal image.  The padding isn't a bug, just a feature
% of working with convolution as you can't get valid values at the edges of
% an image when you convolve, unless you think hard about how to extend
% beyond the image.
theOI = oiCompute(theOI, theScene);
vcAddAndSelectObject(theOI);
oiWindow;

%% Create the coneMosaic object
% Generate default human optics
%defaultHumanOI = oiCreate('wvf human');

% Generate a human foveal cone mosaic.  We construct
% this a bit larger than the image, but at the edges
% you shouldn't trust the computations too much because
% the PSF doesn't know what to convolve with beyond the
% edge of the retinal image.
resamplingFactor = 7;
integrationTimeSecs = 50/1000;
cMosaic = coneMosaicHex(resamplingFactor , ...
    'wave',wls, ...
    'eccBasedConeDensity', false, ...
    'eccBasedConeQuantalEfficiency', false, ...
    'maxGridAdjustmentIterations',200, ...
    'integrationTime', integrationTimeSecs, ...
    'fovDegs', [1.2*backgroundSizeDegs 1.2*backgroundSizeDegs]);

%% Adjust components that determine cone spectral sensitivity
%
% Macular density.  Parallels exposure of lens density in the optical
% image above.  See there for an example of how to adjust this.
% See also "help Macular".
macular0 = cMosaic.macular;
if (any(macular0.wave ~= wls))
    error('Wavelengths in Macular object are not as expected');
end
macularUnitDensity1 = macular0.unitDensity;
macularDensity1 = macular0.density;
macular1 = Macular('wave',wls','unitDensity',macularUnitDensity1,'density',macularDensity1);
cMosaic.macular = macular1;

% Adjust photopigment properties.  Again, similar to how we expose and
% adjust lens and macular pigment density.  See "help receptorPigment"
photopigment0 = cMosaic.pigment;
if (any(photopigment0.wave ~= wls))
    error('Wavelengths in photoPigment object are not as expected');
end
absorbance1 = photopigment0.absorbance;
peakOpticalDensity1 = photopigment0.opticalDensity;
quantalEfficiency1 = photopigment0.peakEfficiency;

% Set this to true set absorbance of photopigments to zero near the test wavelengths,
% which demonstrates that an adjustment has an effect.  We only demonstrate
% this for testWl == 550, at 840 we get an error from the routine that
% tries to scale the image.
changePhotopigment = false;
if (changePhotopigment & testWl == 550)
    wlsChange = find(testSpdRadiance > 0.001*max(testSpdRadiance));
    absorbance1(wlsChange,:) = zeros(size(absorbance1(wlsChange,:)));
end

% Put adjusted photopigment back into mosaic.
photopigment1 = photoPigment('wave',wls,'opticalDensity',peakOpticalDensity1, ...
    'peakEfficiency',quantalEfficiency1,'absorbance',absorbance1);
cMosaic.pigment = photopigment1;

% Eye movement path.  Here, just one time point and no motion.
nTrialsNum = 1;
emPath = zeros(nTrialsNum, 1, 2);

% Compute isomerizations for each eye position.
cMosaic.compute(theOI);
cMosaic.window;