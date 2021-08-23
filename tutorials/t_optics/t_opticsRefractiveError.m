% Examine the effect of different optics parameters on retinal image
%
% Description:
%    Generates an ISETBio scene from a jpg file. Passes the scene via
%    optics in which the user controls several optical factors, such as
%    pupil size, retinal eccentricity, and refractive error.
%

% History:
%    08/24/21       NPC Wrote it. Copyright, ISETBIO Team, 2021

% Your favorite image source
fname = 'ma_blkbkjackal_412.jpg';

% Mean luminance of the scene
meanLuminanceCdPerM2 = 100;

% Field of view of the scene
fieldOfViewDegs = 10;

% Generate an ISETBIo hyperspectral scene assuming it is displayed on a
% typical Apple LCD display

scene = sceneFromFile(fname, 'rgb', meanLuminanceCdPerM2, 'LCD-Apple.mat');
scene = sceneSet(scene, 'wAngular', fieldOfViewDegs);

% Optics at a given eccentricity
eccentricityDegs = [0 0];

% Pupil diameter in MM
pupilDiameterMM = 3.0;

% Human retina magnification factor (retinal microns/deg of visual angle)
retinalMagnificationMicronsPerDegree = 300;

% Additional refractive error in Diopters
additionalRefractiveErrorDiopters = 0.0;

% Generate optics based on human wavefront aberration measurements and the
% subject's normal refractive error
theOI = PolansOptics.oiForSubjectAtEccentricity(1, ...
      'right eye', eccentricityDegs, pupilDiameterMM, ...
      sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
      'refractiveErrorDiopters', additionalRefractiveErrorDiopters);

  
% Compute the retinal image by passing the scene object via the optical image
% object
theOINormalRefractiveError = oiCompute(scene,theOI);

% Generate optics based on human wavefront aberration measurements and a
% an additional refractive error
additionalRefractiveErrorDiopters = 3.0;
theOI = PolansOptics.oiForSubjectAtEccentricity(1, ...
      'right eye', eccentricityDegs, pupilDiameterMM, ...
      sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
      'refractiveErrorDiopters', additionalRefractiveErrorDiopters);
  
theOIAdditionalRefractiveError = oiCompute(scene,theOI);

% Visualize scene and its retinal image under the employed optics
figure()
subplot(3,1,1);
image(sceneGet(scene, 'rgb image'));
axis 'image'
title(sprintf('scene (FOV: %2.2f degrees)', fieldOfViewDegs));

subplot(3,1,2);
image(oiGet(theOINormalRefractiveError, 'rgb image'));
axis 'image'
title(sprintf('retinal image\nnormal refractive error'));

subplot(3,1,3);
image(oiGet(theOIAdditionalRefractiveError, 'rgb image'));
axis 'image'
title(sprintf('retinal image\nadditional refractive error: %2.2f D', additionalRefractiveErrorDiopters));
