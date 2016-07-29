%% t_coneMosaicHex4
%
% Computes hex mosaic isomerization maps for an achromatic Gabor scene and
% illustrates coneMosaicHex's own method for mosaic activation visualization.
%
% NPC, ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

rng('default'); rng(219347);

% Generate a hex mosaic with a medium resamplingFactor
mosaicParams = struct(...
      'resamplingFactor', 8, ...
        'spatialDensity', [0 0.62 0.31 0.07],...
             'noiseFlag', false ...
    );
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);

% Set the mosaic's FOV to a wide aspect ratio
theHexMosaic.setSizeToFOVForHexMosaic([0.9 0.6]);
theHexMosaic.displayInfo();


%% Load an achromatic Gabor scene
[dirName,~] = fileparts(which(mfilename()));
load(fullfile(dirName,'GaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene,'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene,oi);  

% Compute isomerizations
tic
fprintf('\nComputing isomerizations ...');
isomerizations = theHexMosaic.compute(oi,'currentFlag',false);
fprintf('Isomerization computation took %2.1f seconds\n', toc);

% Display isomerizations using coneMosaicHex's own 
% mosaic activation visualization method
tic
fprintf('\nVisualizing responses ... ');
theHexMosaic.visualizeActivationMaps(...
    isomerizations, ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'figureSize', [1550 950], ...
    'mapType', 'modulated hexagons', ...   % choose between 'density plot', 'modulated disks and 'modulated hexagons''
    'colorMap', jet(1024) ...
    )
fprintf('Isomerization visualization took %2.1f seconds\n', toc);
