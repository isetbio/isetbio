%% t_coneMosaicHex0
%
% Shows how to generate a custom hexagonal mosaic, including an S-cone free
% region, and a desired S-cone spacing
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

% Freeze random number generator
%rng('default'); rng(219347);

%% Unit test 1: generate a custom hex mosaic

% Mosaic Parameters
mosaicParams = struct(...
    'resamplingFactor', 10, ...                     % controls the accuracy of the hex mosaic grid
    'spatiallyVaryingConeDensity', true, ...        % whether to have an eccentricity based, spatially - varying density
    'customLambda', [], ...                         % custom spacing?
    'centerInMM', [0.0 0.0], ...                    % mosaic eccentricity
    'size', [96 96], ...                            % generate from a rectangular mosaic of 11 x 16 cones
    'spatialDensity', [0 6/10 3/10 1/10]...         % with a LMS density of of 6:3:1
    );

theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.spatiallyVaryingConeDensity, mosaicParams.customLambda, ...
    'name', 'the hex mosaic', ...
    'center', mosaicParams.centerInMM*1e-3, ...
    'size', mosaicParams.size, ...
    'spatialDensity', mosaicParams.spatialDensity ...
    );
theHexMosaic.visualizeGrid();

theHexMosaic.reassignConeIdentities(...
    'sConeMinDistanceFactor', 3.0, ...   % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 45);       % 45/300 = 0.15, so S-cone free radius of 0.15 deg (0.3 degs diameter)

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid('generateNewFigure', true);

