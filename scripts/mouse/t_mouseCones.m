% t_mouseCones
%
%  Illustrate an image on the mouse cone mosaic.  Just the spatial
%  sampling.  We don't have the pigments or any of that right yet.  This
%  script can get built up to do that.
%
%  We need to do the rods, of course.  We need to check everything.
% 
%
%%
ieInit


%% Make a scene and oi

fov = 100; % deg

nRays = 8; imageSize = 512;
scene = sceneCreate('rings rays',nRays,imageSize);
scene = sceneSet(scene,'fov',fov);
sceneWindow(scene);

%% Not working because of um per degree issue
%{
wvfM  = wvfCreate;
wvfM = wvfSet(wvfM,'um per degree',35);
wvfM = wvfSet(wvfM,'calc pupil size',1);

wvfM  = wvfComputePSF(wvfM);
oiM   = wvf2oi(wvfM);

oiM = oiCompute(oiM,scene);
oiM = oiSet(oiM,'name','mouse');
oiWindow(oiM);

%% Not working correctly
oiGet(oiM,'distance per degree','microns')
oiGet(oiM,'power')

oiPlot(oiM,'psf 550');
%}

%% From the Geng/Dubra/Williams
% The number of microns per degree is about right.
%
% oi = oiCreate('diffraction limited');  % Makes a shift invariant
oi = oiCreate('human');  % Makes a shift invariant

oi = oiSet(oi,'optics fnumber',1);
oi = oiSet(oi,'optics focal length',1.9*1e-3);
oi = oiSet(oi,'lens density',0.0);

% We need to check the PSF and do better here.  Probably we need to deal
% with chromatic aberration, too.
oi = oiCompute(oi,scene);
oiWindow(oi);

%% Check this number as per Geng et al. 34 microns

oiGet(oi,'distance per degree','microns')
oiGet(oi,'power')

oiPlot(oi,'psf 550');

%% Create and visualize a mouse cone mosaic

center = [10 10];
fov = [2 2];
% [~,params] = cMosaic('params');
cmM = cMosaic('eccentricityDegs',center,'sizeDegs',fov);
cmM.visualize;

noiseFree = cmM.compute(oi);

% {
vParams = cm.visualize('params');
vParams.activation = noiseFree;
vParams.activationColorMap = hot(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);
%}