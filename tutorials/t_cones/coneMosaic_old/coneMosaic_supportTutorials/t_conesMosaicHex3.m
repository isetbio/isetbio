%%t_conesMosaicHex3
%
% Generates a rectangular and a hexagonal mosaic and visualizes
% isomerization maps of the two mosaics for three different scenes.
%
% NPC ISETBIO Team, Copyright 2016

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

% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 5 listed here.
%    spatiallyVaryingConeDensity - Boolean. Whether to have an eccentricity
%                                  based, spatially - varying density.
%    centerInMM                  - Vector. (X, Y) Center coordinates in mm.
%    size                        - Vector. Generate from a rectangular
%                                  mosaic of specified rows and columns.
%                                  The listed dimensions are 64 x 64 cones.
%    spatialDensity              - Vector. The spatial density vector, 
%                                  containing the KLMS densities. This is
%                                  specified below as a LMS ratio of
%                                  0.62:0.31:0.07 with no Black included.
%    
mosaicParams = struct('resamplingFactor', 5, ...
    'spatiallyVaryingConeDensity', false, ...
    'centerInMM', [0.0 0.0], ...
    'size', [64 64], ...
    'spatialDensity', [0 0.62 0.31 0.07], ...
    'noiseFlag', 'none');
    
%% Unit Test 1:
% Generate two mosaics, a rectangular one and a hexagonal one. The hex
% mosaic is derived from the rectangular mosaic instantiated below
theRectMosaic = coneMosaic('name', 'the rect mosaic', ...
    'size', mosaicParams.size, ...
    'spatialDensity', mosaicParams.spatialDensity, ...
    'noiseFlag', mosaicParams.noiseFlag);

% Generate a hex mosaic using the pattern of the rectangular mosaic
customLambda = [];
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', customLambda, ...
    'size', mosaicParams.size, ...
    'pattern', theRectMosaic.pattern, ...
    'spatialDensity', mosaicParams.spatialDensity, ...
    'noiseFlag', mosaicParams.noiseFlag, ...
    'name', 'the hex mosaic');
theHexMosaic.displayInfo();

%% Unit Test 2:
% Use the coneMosaic superclass plotting routine in order to visualize the
% two mosaics.
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to display the two mosaics ", ...
        "using the superclass plotting routines. </strong>"));
    pause
end
hFig1 = figure(1);
set(hFig1, 'Position', [10 100 400 400], 'Name', 'Rect mosaic');
hFig2 = figure(2);
set(hFig2, 'Position', [500 100 400 400], 'Name', 'Hex mosaic');
theRectMosaic.plot('cone mosaic', 'hf', hFig1);
theHexMosaic.plot('cone mosaic', 'hf', hFig2);

%% Unit Test 3:
% Compare isomerizations between the rectangular and hexagonal mosaics for
% the ring rays scene.
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to compare isomerizations ", ...
        "between the rectangular and hexagonal mosaics for the ring ", ...
        "rays scene. </strong>"));
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
isomerizationsRect = theRectMosaic.compute(oi, 'currentFlag', false);
isomerizationsHex = theHexMosaic.compute(oi, 'currentFlag', false);
nonNullConeIndices = theHexMosaic.pattern > 1;
allIsomerizations = ...
    [isomerizationsRect(:); isomerizationsHex(nonNullConeIndices)];
isomerizationsRange = prctile(allIsomerizations, [5 95]);


% Render activation images for the hex mosaic
[activationsHexImage, activationsLMSHexImage] = ...
    theHexMosaic.computeActivationDensityMap(isomerizationsHex);


% Display results
hFig = figure(100);
clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101);
clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 1)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 2)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 3)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));

%% Unit Test 4:
% Compare isomerizations between the rectangular and hexagonal mosaics for
% the Vernier scene.
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to compare isomerizations ", ...
        "between the rectangular and hexagonal mosaics for the ", ...
        "Vernier scene. </strong>"));
    pause
end

% Generate the vernier scene
scene = sceneCreate('vernier');
scene.distance = 1;
scene = sceneSet(scene, 'fov', 1.0);
vcAddObject(scene);
sceneWindow

% Compute the optical image
oi = oiCompute(scene, oi);

% Compute isomerizations for both mosaics
isomerizationsRect = theRectMosaic.compute(oi, 'currentFlag', false);
isomerizationsHex = theHexMosaic.compute(oi, 'currentFlag', false);
allIsomerizations = ...
    [isomerizationsRect(:); isomerizationsHex(nonNullConeIndices)];
isomerizationsRange = prctile(allIsomerizations, [5 95]);

% Render activation images for the hex mosaic
[activationsHexImage, activationsLMSHexImage] = ...
    theHexMosaic.computeActivationDensityMap(isomerizationsHex);

% Display results
hFig = figure(100);
clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101);
clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 1)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 2)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 3)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));



%% Unit Test 5:
% Compare isomerizations between the rectangular and hexagonal mosaics for
% the slanted bar scene.
if (interactiveMode)
    commandwindow
    fprintf(strcat("\n<strong>Hit enter to compare isomerizations ", ...
        "between the rectangular and hexagonal mosaics for the ", ...
        "slanted bar scene.</strong>\n\n"));
    pause
end

% Generate the slanted bar scene
scene = sceneCreate('slanted bar');
scene = sceneSet(scene, 'fov', 1.0);
vcAddObject(scene);
sceneWindow;

% Compute the optical image
oi = oiCompute(scene, oi);

% Compute isomerizations for both mosaics
isomerizationsRect = theRectMosaic.compute(oi, 'currentFlag', false);
isomerizationsHex = theHexMosaic.compute(oi, 'currentFlag', false);
allIsomerizations = ...
    [isomerizationsRect(:); isomerizationsHex(nonNullConeIndices)];
isomerizationsRange = prctile(allIsomerizations, [5 95]);

% Render activation images for the hex mosaic
[activationsHexImage, activationsLMSHexImage] = ...
    theHexMosaic.computeActivationDensityMap(isomerizationsHex);

% Display results
hFig = figure(100);
clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101);
clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], ...
    'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 1)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 2)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, ...
    squeeze(activationsLMSHexImage(:, :, 3)));
set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));
