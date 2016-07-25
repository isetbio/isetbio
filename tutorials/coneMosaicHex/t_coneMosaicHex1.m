%% t_coneMosaicHex1
%
% Shows how to change the FOV.
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;
    
rng('default'); rng(219347);

%% Unit test: change FOV
% Mosaic Parameters
mosaicParams = struct(...
      'resamplingFactor', 5, ...
                  'size', [11 16], ...
        'spatialDensity', [0 1/3 1/3 1/3]...
    );
% Generate the hex mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
                   'name', 'the hex mosaic', ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity ...
    );
theHexMosaic.displayInfo();


newFOV = [0.4 0.4];
commandwindow
fprintf('Hit enter to change FOV to [%2.2f, %2.2f]\n', newFOV(1), newFOV(2)); pause
theHexMosaic.setSizeToFOVForHexMosaic(newFOV);
theHexMosaic.displayInfo();

newFOV = [0.3 0.4];
commandwindow
fprintf('Hit enter to change FOV to [%2.2f, %2.2f]\n', newFOV(1), newFOV(2)); pause
theHexMosaic.setSizeToFOVForHexMosaic(newFOV);
