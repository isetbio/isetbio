%% t_coneMosaicHex4
%
% Demonstrates jitter ?
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

rng('default'); rng(219347);

mosaicParams = struct(...
      'resamplingFactor', 4, ...
                  'size', [10 10], ...
        'spatialDensity', [0 0.62 0.31 0.07],...
             'noiseFlag', false ...
    );

% Generate a hex mosaic using the pattern of the Rect mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);
theHexMosaic.setSizeToFOVForHexMosaic([0.5 0.5]);
theHexMosaic.displayInfo();


%% Load an achromatic Gabor scene
[dirName,~] = fileparts(which(mfilename()));
load(fullfile(dirName,'GaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene,'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene,oi);  

% Compute isomerization
isomerizations = theHexMosaic.compute(oi,'currentFlag',false);

% Visualize isomerization maps
theHexMosaic.visualizeActivationMaps(isomerizations);

