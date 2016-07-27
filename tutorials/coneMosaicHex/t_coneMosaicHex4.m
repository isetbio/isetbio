%% t_coneMosaicHex4
%
% Shows hex mosaic isomerization maps for an achromatic Gabor scene 
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

rng('default'); rng(219347);

mosaicParams = struct(...
      'resamplingFactor', 4, ...
        'spatialDensity', [0 0.62 0.31 0.07],...
             'noiseFlag', false ...
    );

% Generate a hex mosaic using the pattern of the Rect mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);
tic
fprintf('\nResising ....');
theHexMosaic.setSizeToFOVForHexMosaic([1 1]);
fprintf('Mosaic resizing took %2.1f seconds\n', toc);
theHexMosaic.displayInfo();


%% Load an achromatic Gabor scene
[dirName,~] = fileparts(which(mfilename()));
load(fullfile(dirName,'GaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene,'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene,oi);  

% Compute isomerization
tic
fprintf('\nComputing isomerizations ...');
isomerizations = theHexMosaic.compute(oi,'currentFlag',false);
fprintf('Isomerization computation took %2.1f seconds\n', toc);

% Visualize isomerization maps
tic
fprintf('\nVisualizing responses ... ');
theHexMosaic.visualizeActivationMaps(isomerizations);
fprintf('Isomerization visualization took %2.1f seconds\n', toc);
