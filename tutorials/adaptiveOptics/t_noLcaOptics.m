% Show how to create an oi structure with no LCA
%
% Description:
%   ISETBio's wavefront optics normally includes wavelength-dependent
%   defocus to model longitudinal chromatic aberration (LCA). This tutorial shows
%   how to construct an optical image object that does not do that.  Useful
%   for modeling experiments where LCA has been corrected.
%
%   The option for turning off LCA is implemented fairly far down the set
%   of structures that manage the optics, and this tutorial shows how to
%   access the right level.
%
%   Retinal images are produced for diffraction limited seeing, with and
%   without LCA.  For the case with, stimuli at the accommodated wavelength
%   are diffraction limited, while at other wavelengths a defocus terms is
%   added to the optics specification according to typical human LCA.
%
%   You can adjust the accommodated wavelength and see the effect of LCA in
%   one of the computed retinal images and not the other.  The test scene
%   used employs wavelengths of 530 and 660.
%
%   The diffraction limited PSF does depend on wavelength, but that effect
%   is small compared to LCA, and is hard to see here.
%
%   You could explore adding higher order aberrations to the
%   wavefront function and looking at how they effect PSF with wavelength.

% History:
%  05/01/19  dhb  Version with LCA defeated.
%  12/xx/20  dhb  Made this self-contained and moved from where it was
%                 tucked away into ISETBio itself.

%% Initialize
close all;

%% Parameters
%
% Make a zero vector of Zernike coefficients to
% represent a diffraction limited pupil function,
% and a few other things.
pupilDiameterMm = 6;
wave = (400:10:700)';
accommodatedWavelength = 800;
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

%% Make a scene with two lines
%
% Read a display
presentationDisplay = displayCreate('CRT12BitDisplay');

% Value to put in each gun as background (on 0-1 scale)
background_level=0.1;

% Image spatial properties, in pixels
imageSize=140;
spacing = 6;
line_thickness=4;
line_height=20;

% Build two line image with separate red and green lines.
%
% First build separate red and green image planes
img_center = floor([imageSize/2,imageSize/2]);
img        = ones(imageSize)*background_level;
img_red    = img;
img_red(img_center - line_height/2:img_center + line_height/2,...
    img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 1;
img_green  = img;
img_green(img_center - line_height/2:img_center + line_height/2,...
    img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 1;
 
% Now put image planes into a fresh image
twoLineImage=zeros(imageSize,imageSize,3);
twoLineImage(:, :, 1) = img_red;
twoLineImage(:, :, 2) = img_green;

% Convert image to ISETBio scene
scene = sceneFromFile(twoLineImage, 'rgb', [], presentationDisplay);
scene = sceneSet(scene, 'fov', 0.5);
visualizeScene(scene);

%% Compute and visualize the retinal images with and without LCA
theOINoLca = oiCompute(theOINoLca, scene);
visualizeOpticalImage(theOINoLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

theOIWithLca = oiCompute(theOIWithLca, scene);
visualizeOpticalImage(theOIWithLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

