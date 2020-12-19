% Show how to create an oi structure with no LCA
%
% Description:
%   ISETBio's wavefront optics normally includes wavelength-dependent
%   defocus to model longitudinal chromatic aberration (LCA). This tutorial shows
%   how to construct an optical image object that does not do that.  Useful
%   for modeling experiments where LCA has been corrected for.
%
%   Retinal images are produced for diffraction limited seeing, with and
%   without LCA.  For the case with, stimuli at the accommodated wavelength
%   are diffraction limited, while at other wavelengths a defocus terms is
%   added to the optics specification according to typical human LCA.
%
%   You can adjust the accommodated wavelength and see the effect of LCA in
%   one of the computed retinal images and not the other.  The test scene
%   used employs wavelengths of 530 and 660.

% History:
%  05/01/19  dhb  Version with LCA defeated.

%% Initialize
close all;

%% Parameters
%
% Make a zero vector of Zernike coefficients to
% represent a diffraction limited pupil function,
% and a few other things.
pupilDiameterMm = 6;
wave = (400:10:700)';
accommodatedWavelength = 530;
zCoeffs = zeros(66,1);

%% Set up wavefront optics object
%
% Compute pupil function using 'no lca' key/value pair to turn off LCA.
% You can turn it back on to compare the effect.
wvfP = wvfCreate('calc wavelengths', wave, 'zcoeffs', zCoeffs, ...
    'name', sprintf('human-%d', pupilDiameterMm));
wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMm);
wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMm);

% Deal with best focus by specifying that the wavefront parameters
% were measured at the wavelength we want to say is in focus. This
% is a little bit of a hack but seems OK for the diffraction limited case
% we're using here.
wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWavelength);

%% Make optical image object using wvfP and no LCA calc
%
% Same as above but don't defeat LCA calc
wvfPNoLca = wvfComputePupilFunction(wvfP,false,'no lca',true);
wvfPNoLca = wvfComputePSF(wvfPNoLca);
theOINoLca = wvf2oi(wvfPNoLca);
opticsNoLca = oiGet(theOINoLca, 'optics');
opticsNoLca = opticsSet(opticsNoLca, 'model', 'shift invariant');
opticsNoLca = opticsSet(opticsNoLca, 'name', 'human-wvf-nolca');
theOINoLca = oiSet(theOINoLca,'optics',opticsNoLca);

%% For comparison, make optical image object using wvfP with LCA calc
wvfPWithLca = wvfComputePupilFunction(wvfP,false,'no lca',false);
wvfPWithLca = wvfComputePSF(wvfPWithLca);
theOIWithLca = wvf2oi(wvfPWithLca);
opticsWithLca = oiGet(theOIWithLca, 'optics');
opticsWithLca = opticsSet(opticsWithLca, 'model', 'shift invariant');
opticsWithLca = opticsSet(opticsWithLca, 'name', 'human-wvf-withlca');
theOIWithLca = oiSet(theOIWithLca,'optics',opticsWithLca);

%% Get a scene.
presentationDisplay = displayCreate('CRT12BitDisplay');
scene = generateTwoLineScene(presentationDisplay, 1, 6); % 1=RG
scene = sceneSet(scene, 'fov', 0.5);
visualizeScene(scene);

%% Compute and visualize the retinal images with and without LCA
theOINoLca = oiCompute(theOINoLca, scene);
visualizeOpticalImage(theOINoLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

theOIWithLca = oiCompute(theOIWithLca, scene);
visualizeOpticalImage(theOIWithLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

% Experimental findings:
%
% These scenes should be distinct in space (clearly 2 lines), but colors are all the same:
% scene = generateTwoLineScene(presentationDisplay,1,20); % RG
% scene = generateTwoLineScene(presentationDisplay,2,20); % GR
% scene = generateTwoLineScene(presentationDisplay,3,20); % YY
%
% Color is only reliable at much greater spacings (>50 pixels)
