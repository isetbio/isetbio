% Hex mosaid isomerization maps for achromatic Gabor scene & Vernier space.
%
% Description:
%    Computes hex mosaic isomerization maps for an achromatic Gabor scene
%    and for the Vernier scene and illustrates coneMosaicHex's own method
%    for mosaic activation visualization.
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

% Freeze rng
rng('default');
rng(219347);

% Generate a hex mosaic with a medium resamplingFactor, using the following
% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 9 listed here.
%    spatiallyVaryingConeDensity - Boolean. Whether to have an eccentricity
%                                  based, spatially - varying density.
%    customLambda                - Numeric. Custom lambda for spacing.
%    spatialDensity              - Vector. The spatial density vector, 
%                                  containing the KLMS densities. This is
%                                  specified below as a LMS ratio of
%                                  0.62:0.31:0.07 with no Black included.
%    noiseFlag                   - String. Type of noise flag.
mosaicParams = struct('resamplingFactor', 9, ...
    'spatiallyVaryingConeDensity', false, ...
    'customLambda', [], ...
    'spatialDensity', [0 0.62 0.31 0.07], ...
    'noiseFlag', 'none');

theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'fovDegs', [0.9 0.6], ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'spatialDensity', mosaicParams.spatialDensity, ...
    'noiseFlag', mosaicParams.noiseFlag);

% Set the mosaic's FOV to a wide aspect ratio
theHexMosaic.displayInfo();

%% Unit Test 1:
% Achromatic Gabor scene
[dirName, ~] = fileparts(which(mfilename()));
load(fullfile(dirName, 't_conesGaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene, 'fov', 1.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene, oi);

% Compute isomerizations
tic
fprintf('\nComputing isomerizations ...');
isomerizationsGabor = theHexMosaic.compute(oi, 'currentFlag', false);
fprintf('Isomerization computation took %2.1f seconds\n', toc);

%% Display isomerizations using mosaic's activation visualization method
tic
fprintf('\nVisualizing responses ... ');
% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap   - Map. The color map to use when displaying the
%                 activation level.
%    figureSize - Vector. The figure size in pixels.
theHexMosaic.visualizeActivationMaps(isomerizationsGabor, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', jet(1024), ...
    'figureSize', [1550 950]);
fprintf('Isomerization visualization took %2.1f seconds\n', toc);

%% Unit Test 2:
% The Vernier scene
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to visualize the hex mosaic ", ...
        "activation maps for the vernier scene. </strong>"));
    pause
end

% Generate the vernier scene
scene = sceneCreate('vernier');
scene.distance = 1;
scene = sceneSet(scene, 'fov', 1.0);
vcAddObject(scene);
sceneWindow

% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

% Compute isomerizations for both mosaics
isomerizationsVernier = theHexMosaic.compute(oi, 'currentFlag', false);

fprintf('\nVisualizing responses ... ');
% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap   - Map. The color map to use when displaying the
%                 activation level.
%    figureSize - Vector. The figure size in pixels.
theHexMosaic.visualizeActivationMaps(isomerizationsVernier, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', bone(1024), ...
    'figureSize', [1550 950]);
fprintf('Isomerization visualization took %2.1f seconds\n', toc);

%% Unit Test 3:
% The rays scene
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to visualize the hex mosaic ", ...
        "activation maps for the rays scene. </strong>"));
    pause
end

% Generate ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene, 'fov', 1.0);
vcAddObject(scene);
sceneWindow

% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

% Compute isomerizations for both mosaics
isomerizationsRays = theHexMosaic.compute(oi, 'currentFlag', false);

fprintf('\nVisualizing responses ... ');
% Activation Map Parameters:
%    Response   - Matrix. The response matrix.
%    mapType    - String. Cone display method, choosing between 'density
%                 plot', 'modulated disks', and 'modulated hexagons'.
%    signalName - String. The name of the signal.
%    colorMap   - Map. The color map to use when displaying the
%                 activation level.
%    figureSize - Vector. The figure size in pixels.
theHexMosaic.visualizeActivationMaps(isomerizationsRays, ...
    'mapType', 'modulated hexagons', ...
    'signalName', 'isomerizations (R*/cone/integration time)', ...
    'colorMap', bone(1024), ...
    'figureSize', [1550 950]);
fprintf('Isomerization visualization took %2.1f seconds\n', toc);
