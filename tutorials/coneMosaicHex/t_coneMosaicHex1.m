%% t_coneMosaicHex1
%
% Shows how to generate a default hexagonal mosaic, and how to customize it (FOV, resamplingFactor).
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

% Interactive mode. Set to true to have it pause at useful places.
% Default is false so we can autopublish without user input
interactiveMode = false;

% Freeze random number generator      
rng('default'); rng(219347);

%% Unit test 1: generate a hex mosaic using defaults params of the superclass (coneMosaic)
% Mosaic Parameters
mosaicParams = struct(...
            'resamplingFactor', 3,  ...        % controls the accuracy of the hex mosaic grid
 'spatiallyVaryingConeDensity', false ...      % whether to have an eccentricity based, spatially - varying density
  );
% Generate the hex grid
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.spatiallyVaryingConeDensity); 
% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();

%% Unit test 2: generate a custom hex mosaic
if (interactiveMode)
    commandwindow;
    fprintf('\n<strong>Hit enter to generate a customized hex mosaic based on an 11x16 rect mosaic with equal LMS proportions\n</strong>');
    fprintf('<strong>Here we use a high resamplingFactor (10) to get a near perfect hex grid\n</strong>');
    pause
end

% Mosaic Parameters
mosaicParams = struct(...
            'resamplingFactor', 10, ...              % controls the accuracy of the hex mosaic grid
'spatiallyVaryingConeDensity', false, ...            % whether to have an eccentricity based, spatially - varying density
                 'centerInMM', [0.5 0.5], ...        % mosaic eccentricity
                       'size', [11 16], ...          % generate from a rectangular mosaic of 11 x 16 cones
             'spatialDensity', [0 1/3 1/3 1/3]...    % with a LMS density of of 0.33:0.33:0.33
    );
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.spatiallyVaryingConeDensity, ...
                   'name', 'the hex mosaic', ...
                 'center', mosaicParams.centerInMM*1e-3, ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity ...
    );
% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();

%% Unit test 3: change the field of view to 0.4 x 0.4 deg
newFOV = [0.4 0.4];
if (interactiveMode)
    commandwindow
    fprintf('\n<strong>Hit enter to change FOV to [%2.2f, %2.2f]\n</strong>', newFOV(1), newFOV(2));
    pause
end
theHexMosaic.setSizeToFOVForHexMosaic(newFOV);
% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();

%% Unit test 3: change the field of view to 0.2 x 0.4 deg
newFOV = [0.2 0.4];
if (interactiveMode)
    commandwindow
    fprintf('\n<strong>Hit enter to change FOV to [%2.2f, %2.2f]\n</strong>', newFOV(1), newFOV(2));
    pause
end
theHexMosaic.setSizeToFOVForHexMosaic(newFOV);
% Print some grid info and visualize it
theHexMosaic.displayInfo();
theHexMosaic.visualizeGrid();