% How to generate a hexagonal mosaic with a spatially-varying cone density.
%
% Description:
%    This tutorial covers how to generate hexagonal mosaic with
%    spatially-varying cone density.
%

% History:
%    xx/xx/16  NPC  ISETBIO Team, Copyright 2016
%    09/14/18  jnm  Formatting

%% Initialize
ieInit;
clear;
close all;

% Interactive mode. Set to true to have it pause at useful places.
% Default is false so we can autopublish without user input
interactiveMode = false;

% Freeze random number generator
rng('default');
rng(219347);

%% Unit Test 1:
% Generate a hex mosaic with spatially-varying cone density positioned at
% (0.5mm, 0.0mm) 

% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 7 listed here.
%    spatiallyVaryingConeDensity - Boolean. Whether to have an eccentricity
%                                  based, spatially - varying density.
%    customLambda                - Numeric. Custom lambda for spacing.
%    size                        - Vector. The mosaic size. Here is 60x60.
%    noiseFlag                   - String. Type of noise flag.
mosaicParams = struct('resamplingFactor', 7, ...
    'spatiallyVaryingConeDensity', true, ...
    'customLambda', [], ...
    'size', [60 60]);
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to create a hex mosaic with ", ...
        "spatially-varying cone density positioned at x = %2.2f mm, ", ...
        "y = %2.2fmm. </strong>"), mosaicParams.centerInMM(1), ...
        mosaicParams.centerInMM(2));
    pause
end

% Generate the hex grid
theHexMosaic1 = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'size', mosaicParams.size, ...
    'maxGridAdjustmentIterations', 150);
% Print some grid info and visualize it
theHexMosaic1.displayInfo();

% Grid Parameters
%    overlayConeDensityContour - String. The contour type. Choose between
%                                'measured', 'theoretical', and 'none'. 
%    coneDensityControuLevels  - Vector. When to draw cone density contours
%    generateNewFigure         - Boolean. Whether to make a new figure.
theHexMosaic1.visualizeGrid('overlayConeDensityContour', 'measured', ...
    'conedensitycontourlevels', [40:50:250] * 1000, ...
    'generateNewFigure', true);

% Grid Parameters
%    overlayConeDensityContour - String. The contour type. Choose between
%                                'measured', 'theoretical', and 'none'. 
%    coneDensityControuLevels  - Vector. When to draw cone density contours
%    generateNewFigure         - Boolean. Whether to make a new figure.
theHexMosaic1.visualizeGrid(...
    'overlayConeDensityContour', 'theoretical', ...
    'conedensitycontourlevels', [40:50:250] * 1000, ...
    'generateNewFigure', true);

%% Unit Test 2:
% Generate a hex mosaic with spatially-uniform cone density positioned
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to create a hex mosaic with ", ...
        "spatially-uniform cone density positioned at x = %2.2f mm, ", ...
        "y = %2.2f mm. </strong>"), mosaicParams.centerInMM(1), ...
        mosaicParams.centerInMM(2));
    pause
end
mosaicParams.eccentricityBasedConeDensity = false;
theHexMosaic2 = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.eccentricityBasedConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'size', mosaicParams.size, ...
    'maxGridAdjustmentIterations', 150);

% Print some grid info and visualize it
theHexMosaic2.displayInfo();
theHexMosaic2.visualizeGrid('generateNewFigure', true);

%% Unit Test 3:
% Compute and display activation maps to a Gabor stimulus
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to compute and visualize ", ...
        "isomerizations maps for the 2 mosaics for an achromatic ", ...
        "Gabor scene. </strong>"));
    pause
end

% Load acrhomatic Gabor scene
[dirName, ~] = fileparts(which(mfilename()));
load(fullfile(dirName, 't_conesGaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene, 'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene, oi);

% Compute isomerizations for the different mosaics and display the results
isomerizationsGabor1 = theHexMosaic1.compute(oi, 'currentFlag', false);
isomerizationsGabor2 = theHexMosaic2.compute(oi, 'currentFlag', false);

activationLUT = jet(1024);
% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap   - Map. The color map to use when displaying the
%                 activation level.
%    figureSize - Vector. The figure size in pixels.
theHexMosaic1.visualizeActivationMaps(isomerizationsGabor1, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', activationLUT, ...
    'figureSize', [1550 950]);

% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap     - Map. The color map to use when displaying the
%                 activation level.
%    figureSize   - Vector. The figure size in pixels.
theHexMosaic2.visualizeActivationMaps(isomerizationsGabor2, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', activationLUT, ...
    'figureSize', [1550 950]);

%% Unit Test 4:
% Compute and display activation maps to the Rays stimulus
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to compute and visualize ", ...
        "isomerizations maps for the 2 mosaics for the rays ", ...
        "scene. </strong>"));
    pause
end
% Generate ring rays stimulus
raysScene = sceneCreate('rings rays');
raysScene = sceneSet(raysScene, 'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(raysScene, oi);

% Compute isomerizations for the different mosaics and display the results
isomerizationsRays1 = theHexMosaic1.compute(oi, 'currentFlag', false);
isomerizationsRays2 = theHexMosaic2.compute(oi, 'currentFlag', false);

activationLUT = bone(1024);
% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap     - Map. The color map to use when displaying the
%                 activation level.
%    figureSize   - Vector. The figure size in pixels.
theHexMosaic1.visualizeActivationMaps(isomerizationsRays1, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', activationLUT, ...
    'figureSize', [1550 950]);

% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap     - Map. The color map to use when displaying the
%                 activation level.
%    figureSize   - Vector. The figure size in pixels.
theHexMosaic2.visualizeActivationMaps(isomerizationsRays2, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', activationLUT, ...
    'figureSize', [1550 950]);
