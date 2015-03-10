% s_oiRayTraceCompute
%
% Show oiCompute's ray tracing version.
%
% Walk through the calculations in oiCompute to see how a scene radiance is
% converted through a lens to an optical image (irradiance) using a ray
% trace example.
%
% NOTES:
%  1) This is broken with an error from opticsRayTrace that says "Ray trace routines not found"
%
% Copyright ImagEval Consultants, LLC, 2011.

%% This is the basic radiance to irradiance code 
scene = sceneCreate('point array',384,32);   % Creates an array of points
scene = sceneSet(scene,'fov',20);  

% These computatons can take a while, so we only do them at a few
% wavelengths for this demonstraton.
scene = sceneInterpolateW(scene,(550:100:650));
% vcAddAndSelectObject(scene); sceneWindow;

%% Create the ray trace optics object
oi = oiCreate('ray trace');

% wave = 550;
% fhmm = 0.5;
% rtPlot(oi,'psf',wave,fhmm);
% psfMovie(oiGet(oi,'optics'));
% vcAddAndSelectObject(oi); oiWindow;

% Notice that the optics model type is set to ray trace
optics = oiGet(oi,'optics');
opticsGet(optics,'model')

%% So the oiCompute will call opticsRayTrace to do the computation
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow;

%% Read the s_opticsRTGridLines.m script to see the details.
%  That script unpakcs the calls inside of oiCompute for the ray trace
%  computation. 

%% End
