% v_photonNoise
%
% Check the routines that generate photon noise for oiPhotonNoise.
%
% Copyright 2013, Imageval, LLC

%% Create a uniform scene with few photons.  

% The small number of photons gives us a chance to see the noise
% distribution on the signal.
scene = sceneCreate('uniform');
scene = sceneSet(scene,'fov',20);  % Pretty big
scene = sceneAdjustLuminance(scene,10^-11);
% vcAddAndSelectObject(scene); sceneWindow

%% Create and crop out center of OI
oi = oiCreate; 

% No lens shading
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'cos4th','off');
oi = oiSet(oi,'optics',optics);

oi = oiCompute(oi, scene);
% Reset the rect if you adjust any sizes
% vcAddAndSelectObject(oi); oiWindow
% [oi,rect] = oiCrop(oi);

rect = [6 6 29 29];  % Middle part
oi = oiCrop(oi,rect);
% vcAddAndSelectObject(oi); oiWindow

%% Compare the variance of the photon noise and the mean level
photons = oiGet(oi,'photons');
pMean = photons(:,:,10);  mean(pMean(:));

noisyPhotons = oiGet(oi,'photons with noise');
pNoise = noisyPhotons(:,:,10) - photons(:,:,10);
% vcNewGraphWin; hist(pNoise(:))

% This should be close to 1
t = var(pNoise(:))/mean(pMean(:));
sprintf('Should be near 1:  %f\n',t)
% assert(abs(t - 1) < 0.1)

