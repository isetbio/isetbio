%% t_coneMosaicHex4
%
% Demonstrates jitter ?
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

rng('default'); rng(219347);

mosaicParams = struct(...
      'resamplingFactor', 5, ...
                  'size', [48 48], ...
        'spatialDensity', [0 0.62 0.31 0.07],...
             'noiseFlag', false ...
    );

% Generate a hex mosaic using the pattern of the Rect mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);
theHexMosaic.displayInfo();


%% Unit test:  the ring rays scene
commandwindow
fprintf('\n<strong>Hit enter to compare isomerizations between the rect and hex mosaics for the ring rays scene. </strong>');
pause

% Generate ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov', 1.0);
vcAddObject(scene); sceneWindow
    
% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene,oi);  

% Compute isomerization
isomerizations = theHexMosaic.compute(oi,'currentFlag',false);

% Visualize isomerization maps
theHexMosaic.visualizeActivationMaps(isomerizations);

