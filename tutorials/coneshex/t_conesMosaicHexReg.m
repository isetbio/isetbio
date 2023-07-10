% Create/use hex. mosaic w/ custom fixed cone spacing & inner segment diam.
%
% Description:
%    Generate and use a hexagonal mosaic with fixed (but custom) cone
%    spacing and inner segment diameter.
%
%    This shows how to generate a custom regular hexagonal mosaic,
%    including a desired cone spacing and inner segment diameter. Here we
%    allow S-cones even at the foveola, and a 4-fold spacing between L/M
%    cones. Then show how to compute the isomerizations for this mosaic to
%    a simple stimulus.
%
% See Also:
%   t_coneMosaicHex
%

% History:
%    xx/xx/16  NPC  ISETBIO Team, Copyright 2016
%    08/08/17  NPC  Fixed and cleaned up for updated @coneMosaicHex class
%    07/23/18  JNM  Formatting

%% Initialize
ieInit;
clear;
close all;

%% Set mosaic parameters
% The various mosaic parameters, and their descriptions.
%
%    'name'             - String. The name of the mosaic.
%    'resamplingFactor' - Numeric. Sets underlying pixel spacing; controls
%                         the accuracy of the hex mosaic grid (9 is pretty
%                         good, but slow)
%    'fovDegs'          - Numeric. The FOV in degrees. 0.35 degrees.
%    'customLambda'     - Numeric. The cone spacing in microns.
%    'customInnerSegmentDiameter'
%                       - Numeric. The inner segment diameter in microns.
%    'sConeFreeRadiusMicrons'
%                       - Numeric. Empty meaning no S-cone free region.
%    'eccBasedConeDensity'
%                       - Boolean. Whether to have an eccentricity based,
%                         spatially - varying density.
%    'eccBasedConeQuantalEfficiency'
%                       - Boolean. Whether to have an eccentricity based,
%                         variation in quantal efficiency due to changes in
%                         inner segment aperture and outer-segment length.
%    'eccBasedMacularPigment'
%                       - Boolean. Whether to have an eccentricity based,
%                         variation in the optical density of the macular
%                         pigment
%    'sConeMinDistanceFactor'
%                       - Numeric. Min distance between neighboring S-cones
%                         following the formula f * local cone separation,
%                         which is then used to make the S-cone lattice
%                         semi-regular.
%    'spatialDensity'   - Vector. The vector containing the KLMS spatial
%                         densities, with a LMS density of of 6:3:1.
mosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...
    'fovDegs', 0.35, ...
    'customLambda', 2.0, ...
    'customInnerSegmentDiameter', 1.6, ...
    'sConeFreeRadiusMicrons', [], ...
    'eccBasedConeDensity', false, ...
    'sConeMinDistanceFactor', 4.0, ...
    'spatialDensity', [0, 6 / 10, 3 / 10, 1 / 10]);

%% Generate the mosaic. This takes a little while.
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'name', mosaicParams.name, ...
    'fovDegs', mosaicParams.fovDegs, ...
    'customLambda', mosaicParams.customLambda, ...
    'customInnerSegmentDiameter', mosaicParams.customInnerSegmentDiameter, ...
    'eccBasedConeDensity', mosaicParams.eccBasedConeDensity, ...
    'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...
    'spatialDensity', mosaicParams.spatialDensity);

%% Print some grid info
theHexMosaic.displayInfo();

%% Visualize the mosaic with both the inner segment and geometric area
% Shows the light collecting area (inner segment) and the geometric area

% Choose aperture from: 'both', 'lightCollectingArea', 'geometricArea'
visualizedAperture = 'lightCollectingArea';
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', visualizedAperture, ...
    'generateNewFigure', true);

%% Compute isomerizations to a simple stimulus for the resulting mosaic.
% Generate ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene, 'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

% Compute isomerizations for both mosaics and look at them in a window.
isomerizationsHex = theHexMosaic.compute(oi, 'currentFlag', false);
theHexMosaic.window;

% Get the numerical values of the isomerizations. Pattern values of 1 on
% the underlying grid indicate a location with no cones. Values of 2, 3, 4
% indicate L, M and S cones respectively.
nonNullConeIndices = theHexMosaic.pattern > 1;
allIsomerizations = isomerizationsHex(nonNullConeIndices);
isomerizationsRange = prctile(allIsomerizations, [5, 95]);
