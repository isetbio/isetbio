%%t_conesMosaicHex  Generate and use a hexagonal cone mosaic.
%
% Description:
%   Shows how to generate a custom hexagonal mosaic, including an S-cone free
%   region, and a desired S-cone spacing.
%
%   Then shows how to compute isomerizations for this mosaic to a simple stimulus.
%
% See also: advancedTutorials/t_conesMosaicHex1, ..., advancedTutorials/t_conesMosaicHex6

% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

%% Set mosaic parameters
mosaicParams = struct(...
    'resamplingFactor', 5, ...                      % Sets underlying pixel spacing; controls the accuracy of the hex mosaic grid
    'spatiallyVaryingConeDensity', true, ...        % Whether to have an eccentricity based, spatially - varying density
    'customLambda', [], ...                         % Custom spacing?
    'centerInMM', [0.0 0.0], ...                    % Mosaic eccentricity
    'size', [64 64], ...                            % Generate from a rectangular mosaic of 96 by 96 cones
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

%% Generate the mosaic.  This takes a little while.
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.spatiallyVaryingConeDensity, mosaicParams.customLambda, ...
    'name', 'the hex mosaic', ...
    'center', mosaicParams.centerInMM*1e-3, ...
    'size', mosaicParams.size, ...
    'spatialDensity', mosaicParams.spatialDensity ...
    );

%% Show the mosaic in a figure window.
theHexMosaic.visualizeGrid();

%% Make the S cone submosaic more regular and add an S-cone free region in the fovea
theHexMosaic.reassignConeIdentities(...
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 45);                  % 45 microns.  With 300 microns/degree -> 45/300 = 0.15, so S-cone free radius of 0.15 deg (0.3 degs diameter)

%% Print some grid info and visualize the second mosaic
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid('generateNewFigure', true);

%% Compute isomerizations to a simple stimulus for the resulting mosaic.

% Generate ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov', 1.0);
    
% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene,oi);  

% Compute isomerizations for both mosaics and look at them in a window.
isomerizationsHex = theHexMosaic.compute(oi,'currentFlag',false);
theHexMosaic.window;

% Get the numerical values of the isomerizations. Pattern values of 1 on
% the underlying grid indicate a location with no cones.  Values of 2, 3, 4
% indicate L, M and S cones respectively.
nonNullConeIndices = theHexMosaic.pattern > 1;
allIsomerizations = isomerizationsHex(nonNullConeIndices);
isomerizationsRange = prctile(allIsomerizations, [5 95]);

