% Resolution in the perifovea region
%

%%
ieInit;
close all

%% Rings and rays
%
nRays    = 10;
nSamples = 256;
fov      = 3;   % Total cone fov will be is fov/2 + fov

%% Make three copies of rings and rays

scene0 = sceneCreate('ringsrays', nRays, nSamples);
scene0 = sceneSet(scene0, 'fov', fov);
scene00 = sceneCombine(scene0,scene0);
scene = sceneCombine(scene00, scene0,'direction','horizontal');
sceneWindow(scene);
drawnow;

%% Onward to oi and coneMosaic
%
oi = oiCompute(scene, oi);
% oiWindow(oi);
drawnow;

%% Generate mosaic centered at target eccentricity
hfov = sceneGet(scene,'hfov');
vfov = sceneGet(scene,'vfov');
cm = cMosaic(...
    'sizeDegs', [hfov vfov], ...   % (x,y)
    'eccentricityDegs', [0 0] ... 
    );

% Compute the noise-free excitation response
noiseFreeExcitationResponse = cm.compute(oi);

%% Have a look

cm.visualize( ...
    'domain', 'degrees', ...
    'activation', noiseFreeExcitationResponse, ...
    'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));

%% END
