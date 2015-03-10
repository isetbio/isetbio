%% s_sceneBitdepth
%
% We illustrate using a bit depth of 32 bits
%
% (c) Imageval Consulting, LLC 2012

%%
s_initISET

%%
fov = 2;
scene = sceneCreate('default', 'bit depth', 32);
scene = sceneSet(scene, 'h fov', fov);
scene.data

vcAddAndSelectObject(scene);
sceneWindow;

%%
oi = oiCreate;
oi = oiSet(oi,'bit depth', 32);
oi = oiCompute(oi,scene);
oi.data

vcAddAndSelectObject(oi);
oiWindow;

%% 
sensor = sensorCreate;
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor);
sensorWindow;
%% End