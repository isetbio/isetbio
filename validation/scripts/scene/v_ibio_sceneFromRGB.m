% v_sceneFromRGB_debug0
% Test how isetbio creates a scene from rgb data using sceneFromFile
%
% Function sceneFromFile uses 1) the rgb data and the spectral power
% distributions (spds) of the display color primaries to calculate the
% spectral radiance of the displayed image, and 2) the spectral power
% distribution of the display white point as an estimate of the scene
% illuminant.
%
% Notes about the method
%     If we have display calibration data, we can accurately predict the
%     radiance emitted when an image is rendered on the display. But we
%     need a scene illuminant to estimate scene reflectances. We use the
%     spectral power of the display whitepoint (max r, max g, max b) as an
%     estimate of the scene illuminant. We then calculate reflectances of
%     surfaces in the scene by dividing the scene radiance by the
%     illuminant spd. The surface reflectances will not be accurate, but
%     they will be feasible. And, more importantly, calculating scene
%     reflectances makes it possible to render the scene under a different
%     illuminant.
%
% In this exammple, the scene illuminant is also changed.
%
% For this example, the gamma function of the display is not used.
%
% Copyright ImagEval, 2011

%% Initialize validation
ieInit

% May need to adjust this some day.  But for now, I am putting the
% validation data for ISETCam tools in isetcam/data/validation.
validFile = fullfile(isetRootPath,'data','validation','scene','v_data_sceneFromRGB');
load(validFile,'validationData');
vData = validationData;
vData.whiteXYZ = [108.6659  118.1495  121.0519];

%% Load display calibration data

displayCalFile = 'LCD-Apple.mat';
load(displayCalFile,'d'); dsp = d;
wave = displayGet(dsp,'wave');
spd = displayGet(dsp,'spd');

ieNewGraphWin; plot(wave,spd);
xlabel('Wave (nm)'); ylabel('Energy'); grid on
title('Spectral Power Distribution of Display Color Primaries');

assert(isequal(wave,vData.wave));
assert(isequal(spd,vData.spd));

%% Analyze the display properties: Chromaticity

d = displayCreate(displayCalFile);
whtSPD = displayGet(d,'white spd');
wave   = displayGet(d,'wave');
whiteXYZ = ieXYZFromEnergy(whtSPD',wave);
fig = chromaticityPlot(chromaticity(whiteXYZ));

assert(max(abs(whiteXYZ - vData.whiteXYZ)) < 1e-4);  % Just rounding
assert(isequal(wave,vData.wave1));
assert(isequal(spd,vData.spd1));

%% Read in an rgb file and create calibrated display values
rgbFile = fullfile(isetRootPath,'data','images','rgb','eagle.jpg');
scene = sceneFromFile(rgbFile,'rgb',[],displayCalFile);
sceneShowImage(scene);

assert(isequal(scene.distance,vData.scene.distance));
assert(isequal(scene.wAngular,vData.scene.wAngular));
assert(isequal(scene.spectrum,vData.scene.spectrum));

%% Internal validation
%
% We are having cross-platform issues with the horizontal angular
% extent of the scene varying.  So we'll compute it from first
% principles here.
assert(isequal(sceneGet(scene,'hfov'),sceneGet(vData.scene,'hfov')));
assert(isequal(sceneGet(scene,'size'),sceneGet(vData.scene,'size')));
assert(abs(sceneGet(scene,'mean luminance') - sceneGet(vData.scene,'mean luminance')) < 1e-3);

%{
testByHand.hFovDegreesFromScene = sceneGet(scene,'hfov');
testByHand.displayDPI = displayGet(d,'dpi');
testByHand.displayMetersPerDot = displayGet(d,'meters per dot');
testByHand.distance = sceneGet(scene,'distance');
testByHand.degreesPerDot = displayGet(d,'deg per dot');
testByHand.degreesPerDotCheck = 2*atand(testByHand.displayMetersPerDot/(2*testByHand.distance));
testByHand.hSceneDots = sceneGet(scene,'cols');
testByHand.hFovMeters = testByHand.displayMetersPerDot*testByHand.hSceneDots;
testByHand.hFovDegreeesCheck = testByHand.hSceneDots*testByHand.degreesPerDot;
%}

%% Change the illuminant to 4000 K
bb = blackbody(sceneGet(scene,'wave'),4000,'energy');
scene = sceneAdjustIlluminant(scene,bb);
sceneShowImage(scene);
assert(isequal(bb,vData.bb));

%% END