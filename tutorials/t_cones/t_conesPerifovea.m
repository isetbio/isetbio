% Resolution in the perifovea region
%

%%
ieInit;
close all

%% Rings and rays
%
nRays    = 10;
nSamples = 256;
fov      = 2;

%% Make three copies of rings and rays

scene = sceneCreate('ringsrays', nRays, nSamples);
scene = sceneSet(scene, 'fov', fov);
scene = sceneCombine(scene, scene,'direction','centered');
sceneWindow(scene);

%% Onward to oi and coneMosaic
%
oi = oiCompute(scene, oi);
% oiWindow(oi);

%% Generate mosaic centered at target eccentricity
hfov = sceneGet(scene,'hfov');
vfov = sceneGet(scene,'vfov');
fprintf('Creating cone mosaic %d by %d\n',hfov,vfov);
tic
cm = cMosaic(...
    'sizeDegs', [hfov vfov], ...   % (x,y)
    'eccentricityDegs', [0 0] ... 
    );
toc

% Compute the noise-free excitation response
noiseFreeExcitationResponse = cm.compute(oi);

%% Have a look

cm.visualize( ...
    'domain', 'degrees', ...
    'activation', noiseFreeExcitationResponse, ...
    'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));

%% END
