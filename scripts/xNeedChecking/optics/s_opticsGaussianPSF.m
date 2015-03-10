%s_DemoGaussianPSF.m

% Script to demonstrate the use of a space-invariant Gaussian PSF
%
% This script creates a scene with a point array.
%
% The scene is transformed using diffraction limited methods
% (shift-invariant).  The diffraction limited methods use the f# and focal
% length.
%
% Copyright ImagEval Consultants, LLC, 2006

clear all

wave = (450:100:650);
nWaves = length(wave);

% Create scene
scene = sceneCreate('pointArray',128,32);
scene = sceneInterpolateW(scene,wave);
scene = sceneSet(scene,'hfov',1);
scene = sceneSet(scene,'name','psfPointArray');
vcAddAndSelectObject('scene',scene); sceneWindow;

% Create optical image
%
oi = oiCreate;
oi = oiSet(oi,'spectrum',sceneGet(scene,'spectrum'));

% Calculate Gaussian PSF
%
psfType = 'gaussian';
xyRatio = 3*ones(1,nWaves);
waveSpread = wave/wave(1);
optics = siSynthetic(psfType,oi,waveSpread,xyRatio);
% optics = opticsSet(optics,'model','shiftinvariant');

% Put optics back into oi and display "blurred" optical image
%
oi = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
vcAddAndSelectObject('oi',oi); 
oiWindow;
