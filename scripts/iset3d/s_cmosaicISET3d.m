%% Illustrating how to compute with iset3d and isetbio
%


%%
ieInit;
if ~piDockerExists, piDockerConfig; end

% Read the recipe
thisSE = sceneEye('chess set scaled','human eye','navarro');
%%

thisSE.set('fov',2);
thisSE.set('rays per pixel',32);  % Pretty quick, but not high quality

oiLeft = thisSE.render;  % Render radiance and depth, and then show
oiWindow(oiLeft);

%% Cone mosaic

cm = cMosaic('sizeDegs',[3,3],'eccentricityDegs',[2,0]);
cm.visualize();

%% Compute 

excitations = cm.compute(oiLeft);
params = cm.visualize('params');
