%% s_opticsRTSynthetic
%
% Illustrate a ray trace calculation based on a synthetic ray trace data
% set built by rtSynthetic.
%
% The synthetic point spreads are bivariate normals that increase from
% center to periphery.  Over time, we will add parameters to the function
% that control the bivariate normal growth with field height as well as
% other parameters. 
%
% (c) Imageval Consulting LLC, 2012

%% Create a test scene
scene = sceneCreate('point array',256);
% scene = sceneCreate('sweepfrequency',256);
scene = sceneSet(scene,'h fov',4);
scene = sceneInterpolateW(scene,550:100:650);
vcAddAndSelectObject(scene);
sceneWindow;

%% Build the optical image
oi = oiCreate;
rtOptics = []; spreadLimits = [1 5]; xyRatio = 1.6;
rtOptics = rtSynthetic(oi,rtOptics,spreadLimits,xyRatio);
oi = oiSet(oi,'optics',rtOptics); 
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name','Synthetic-RT-Increasing-Gaussian');
vcAddAndSelectObject(oi); oiWindow;

%% Show the PSF
rtPlot(oi,'psf',550,0);
rtPlot(oi,'psf',550,1);

