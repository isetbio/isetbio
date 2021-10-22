%% 
ieInit;
if ~piDockerExists, piDockerConfig; end

%%
mccR = piRecipeDefault('scene name','MacBethChecker');
mccScene = piWRS(mccR);
mccScene = sceneSet(mccScene,'fov',1);

oi = oiCreate;
oi = oiCompute(oi,mccScene);
oiWindow(oi);

%% Cone mosaic

% Starts the parallel pool, which is a bit annoying IMHO.
cm = cMosaic('sizeDegs',[1,1],'eccentricityDegs',0);
cm.visualize;

activations = cm.compute(oi);

params = cm.visualize('help');

% Just the mosaic
params.activation = [];
cm.visualize(params);

params.activation = activations;
params.activationColorMap = gray(256);
params.verticalActivationColorBar = true;
params.activationRange = [0 max(noisyExcitations)];
params.backgroundColor = [0.3 0.5 0.7];  % Not sure this works
params.visualizedConeAperture = 'lightCollectingArea';
params.domain = 'degrees';

cm.visualize(params);

%% Rotate the scene

mccScene2 = sceneRotate(mccScene,'ccw');
oi2 = oiCompute(oi,mccScene2);
oiWindow(oi2);

params.activation = cm.compute(oi2);
cm.visualize(params);

%%
cm = cMosaic('sizeDegs',[1,1],'eccentricityDegs',[3 0]);
cm.visualize;
params.activation = cm.compute(oi2);
cm.visualize(params);


