% t_FixationLetterE
% 
% This tutorial is a first-pass simplification of the experiments desribed
% in Anderson & Olshausen, Cosyne 2016, and Ratnam, Harmening and Roorda,
% VSS 2015. The purpose is to demonstrate that Isetbio can be used to
% implement an approximation of the stimulus and generate simulated retinal
% ganglion cell spiking responses that could be used for input into the
% model of Anderson & Olshausen.
% 
% Isetbio is a framework for simulating the physiological optics of the eye
% and the dynamis of cone photoreceptor and retinal ganglion cell
% resposnes. The RGC model is based on the paper by Pillow, Chichilnisky,
% et al., Nature, 2008. At this point, it simply operates on the raw image
% data, but it will soon be connected to the output of the photoreceptor
% responses.
% 
% Even though the RGC simulation is based on data collected in primate
% retina at high eccentricity, the parameters of the model can be
% extrapolated to low eccentricity. The RGC model was built for photopic
% light conditions, like the conditions used for the experiments found in
% Ratnam, Harmening & Roorda. 
% 
% This is an easy way to generate an accurate prediction of spiking output
% for a number of retinal ganglion cell mosaics. Each mosaic can be
% constructed to have the properties of a certain type of cell, like on
% parasol, off parasol, on midget, etc. The simulation loses some accuracy
% due to the absence of data at low eccentricities, but it is a start.
% 
% Outline: 
% 1. Build the letter E stimulus movie, with badly "simulated" eye
% movements.
% 2. Build the representation of the raw RGB data with osIdentity (a
% placeholder for when we don't actually simulate the outer segment current
% response).
% 3. Compute the RGC mosaic response and plot the rasters.
% 
% 3/2016 JRG (c) isetbio team

%% Initialize
clear;
ieInit;
 
%% Load image

% Load the fixation letter E stimulus

% Number of frames
params.nSteps = 10;
% Size of image
params.row = 64;
params.col = 64;
% Field of view
params.fov = 0.7;

% Build stimulus
fixationalE = ieFixationalLetterE(params);

%% Show movie of stimulus
figure;
for frame1 = 1:size(fixationalE.sceneRGB,3)
    imagesc(squeeze(fixationalE.sceneRGB(:,:,frame1,:)));
    colormap gray; drawnow;
end
close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('identity');

sceneSize = sceneGet(fixationalE.scene,'size');
retinalPatchWidth = sensorGet(fixationalE.absorptions,'width','m');
% retinalPatchHeight = sensorGet(movingBar.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

timeStep = sensorGet(fixationalE.absorptions,'time interval','sec');
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', fixationalE.sceneRGB);

%% Build the inner retina object

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(os, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','glm','type','on midget');
% innerRetina0.mosaicCreate('model','glm','type','off midget');
% innerRetina0.mosaicCreate('model','glm','type','on parasol');
% innerRetina0.mosaicCreate('model','glm','type','off parasol');
irPlot(innerRetina0,'mosaic');

%% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, os);

innerRetinaRaster = mosaicGet(innerRetina0.mosaic{1},'responseSpikes');
% innerRetinaRaster is N neurons X M neurons X T trials

% irPlot(innerRetina0, 'psth');
% irPlot(innerRetina0, 'linear');
irPlot(innerRetina0, 'raster');