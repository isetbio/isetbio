% Show how to generate and customize a default hex mosaic.
%
% Description:
%    Shows how to generate a default hexagonal mosaic, and how to customize
%    it (FOV, resamplingFactor).
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
% Generate a hexagonal mosaic using the defaults parameters from the
% coneMosaic superclass.

% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 5 listed here.
%    spatiallyVaryingConeDensity - Boolean. Whether to have an eccentricity
%                                  based, spatially - varying density.
%    customLambda                - Numeric. Custom lambda for spacing.
mosaicParams = struct('resamplingFactor', 5, ...
    'spatiallyVaryingConeDensity', false, ...
    'customLambda', []);

% Generate the hex grid
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda);

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', 'lightCollectingArea');
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', 'geometricArea', 'generateNewFigure', true);

%% Unit Test 2a:
% Generate a custom hex mosaic
if (interactiveMode)
    commandwindow;
    fprintf(strcat("\n<strong>Hit enter to generate a customized hex ", ...
        "mosaic based on an 16x16 rect mosaic with custom lambda (3 ", ...
        "microns)\n</strong>"));
    fprintf(strcat("<strong>Here we use a high resamplingFactor (10) ", ...
        "to get a near perfect hex grid\n</strong>"));
    pause
end

% Generate the hex grid
mosaicParams.resamplingFactor = 10;
mosaicParams.customLambda = 3;  % 3 microns
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'size', [16 16]);

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', 'lightCollectingArea');
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', 'geometricArea', 'generateNewFigure', true);

%% Unit test 2b:
% Generate a custom hex mosaic
if (interactiveMode)
    commandwindow;
    fprintf(strcat("\n<strong>Hit enter to generate a customized hex ", ...
        "mosaic based on an 16x16 rect mosaic with custom lambda ", ...
        "(6 microns)\n</strong>"));
    fprintf(strcat("<strong>Here we use a high resamplingFactor (10) ", ...
        "to get a near perfect hex grid\n</strong>"));
    pause
end

% Generate the hex grid
customLambda = 6;  % 6 microns
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'size', [16 16]);

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();

%% Unit Test 3:
% Generate a custom hex mosaic
if (interactiveMode)
    commandwindow;
    fprintf(strcat("\n<strong>Hit enter to generate a customized hex ", ...
        "mosaic based on an 11x16 rect mosaic with equal LMS ", ...
        "proportions\n</strong>"));
    fprintf(strcat("<strong>Here we use a high resamplingFactor (10) ", ...
    "to get a near perfect hex grid\n</strong>"));
    pause
end

% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 10 listed here.
%    spatiallyVaryingConeDensity - Boolean. Whether to have an eccentricity
%                                  based, spatially - varying density.
%    customLambda                - Numeric. Custom lambda for spacing.
%    size                        - Vector. Generate from a rectangular
%                                  mosaic of specified rows and columns.
%                                  The listed dimensions are 11 x 16 cones.
%    spatialDensity              - Vector. The spatial density vector,
%                                  containing the KLMS densities. This is
%                                  specified below as a LMS ratio of
%                                  0.33:0.33:0.33 with no Black included.
mosaicParams = struct('resamplingFactor', 10, ...
    'spatiallyVaryingConeDensity', false, ...
    'customLambda', [], ...
    'size', [11 16], ...
    'spatialDensity', [0 1/3 1/3 1/3]);
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda', mosaicParams.customLambda, ...
    'name', 'the hex mosaic', ...
    'size', mosaicParams.size, ...
    'spatialDensity', mosaicParams.spatialDensity);

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();
