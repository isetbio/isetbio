% Visualize the effects of different resampling factors.
%
% Description:
%    Visualizes the effects of different resampling factors.
%       * A resampling factor value of 1 gives a rectangular grid.
%       * A resampling factor of 2 gives pretty bad artifacts with
%         largevinhomogenities.
%       * As the resampling factor is increased beyond 2, inhomogeneities
%         start to disappear and the generated grid approximates a perfect
%         hex grid.
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
% display different aspects of a hex mosaic in a 4 panel display, including
% the hex mosaic, the originating rectangular mosaic, the hex grid used for
% sampling, and the null cones. 

% Mosaic Parameters
%    resamplingFactor            - Numeric. Controls accuracy of the hex
%                                  mosaic grid. Factor of 7 listed here.
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
mosaicParams = struct('resamplingFactor', 7, ...
    'spatiallyVaryingConeDensity', false, ...
    'customLambda', [], ...
    'size', [11 16], ...
    'spatialDensity', [0 1/3 1/3 1/3]);

% Generate the hex mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'eccBasedConeDensity', mosaicParams.spatiallyVaryingConeDensity, ...
    'customLambda',  mosaicParams.customLambda, ...
    'name', 'the hex mosaic', ...
    'size', mosaicParams.size, ...
    'spatialDensity', mosaicParams.spatialDensity);

% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();

% Display the corresponding rectangular grid instead
theHexMosaic.visualizeGrid('panelPosition', [1 2], ...
    'showCorrespondingRectangularMosaicInstead', true);

% Overlay the null sensors as well
theHexMosaic.visualizeGrid('panelPosition', [2 1], ...
    'overlayNullSensors', true);

% Overlap the perfect hexagonal mesh
theHexMosaic.visualizeGrid('panelPosition', [2 2], ...
    'overlayHexMesh', true);

keepGoing = interactiveMode;
while (keepGoing)
    commandwindow
    resamplingFactor = input(sprintf(strcat("\n<strong>Enter a new ", ...
        "resampling factor [>= 1]. A negative exits the loop. ", ...
        "New resampling Factor: </strong>")));
    if (isempty(resamplingFactor)), resamplingFactor = 1; end
    
    if (resamplingFactor <= 0), keepGoing = false; continue; end
    theHexMosaic.resampleGrid(resamplingFactor);
    
    % Print some grid info and visualize it
    theHexMosaic.displayInfo();
    theHexMosaic.visualizeGrid();
    
    % Display the corresponding rectangular grid instead
    theHexMosaic.visualizeGrid('panelPosition', [1 2], ...
        'showCorrespondingRectangularMosaicInstead', true);
    
    % Overlay the null sensors as well
    theHexMosaic.visualizeGrid('panelPosition', [2 1], ...
        'overlayNullSensors', true);
    
    % Overlap the perfect hexagonal mesh
    theHexMosaic.visualizeGrid('panelPosition', [2 2], ...
        'overlayPerfectHexMesh', true);
end
