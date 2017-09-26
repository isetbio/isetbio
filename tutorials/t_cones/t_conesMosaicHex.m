%%t_conesMosaicHex  Generate and use a hexagonal cone mosaic with eccentricity-based cone spacing.
%
% Description:
%   Shows how to generate a custom hexagonal mosaic, with an eccentricity based cone spacing 
%   including an S-cone free region, and a desired S-cone spacing.
%
%   Then shows how to compute isomerizations for this mosaic to a simple stimulus.
%
% See also: advancedTutorials/t_conesMosaicHex1, ..., advancedTutorials/t_conesMosaicHex6

% NPC ISETBIO Team, Copyright 2016
%
% 08/08/17  npc   Fixed and cleaned up for updated @coneMosaicHex class
% 08/17/17  dhb   Changed parameters to make this go faster.
% 09/25/17  npc   Updated to show cone density contours.

%% Initialize
ieInit; clear; close all;

%% Set mosaic parameters
mosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...                      % Sets underlying pixel spacing; controls the rectangular sampling of the hex mosaic grid
    'fovDegs', 0.3, ...                             % FOV in degrees
    'eccBasedConeDensity', true, ...                % Whether to have an eccentricity based, spatially - varying density
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 0.15*300, ...         % Radius of S-cone free retina, in microns (300 microns/deg).
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

% 360 seconds
lowQuality.tolerance1 = 0.01*10;
lowQuality.tolerance2 = 0.001*10;
lowQuality.marginF = 1.3;
lowQuality.mosaicFileName = 'lowQMosaic.mat';
lowQuality.saveMosaic = true;

% 1040 seconds
medQuality.tolerance1 = 0.01*3;
medQuality.tolerance2 = 0.001*3;
medQuality.marginF = 1.5;
medQuality.mosaicFileName = 'medQMosaic.mat';
medQuality.saveMosaic = true;

highQuality.tolerance1 = 0.01;
highQuality.tolerance2 = 0.001;
highQuality.marginF = 1.5;
highQuality.mosaicFileName = 'highQMosaic.mat';
highQuality.saveMosaic = true;

highQuality2.tolerance1 = 0.01;
highQuality2.tolerance2 = 0.001;
highQuality2.marginF = 2.0;
highQuality2.mosaicFileName = 'highQ2Mosaic.mat';
highQuality2.saveMosaic = true;


% Choose between low and high quality params
qParams = lowQuality;
qParams = highQuality2;

tic
%% Generate the mosaic.  This takes a little while.
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'name', mosaicParams.name, ...
    'fovDegs', mosaicParams.fovDegs, ...
    'eccBasedConeDensity', mosaicParams.eccBasedConeDensity, ...
    'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ... 
    'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...                   
    'spatialDensity', mosaicParams.spatialDensity, ...
    'latticeAdjustmentPositionalToleranceF', qParams.tolerance1, ...   % This value is too high and is chosen to reduce the compute time. For production work, either do not set, or set to equal or lower than 0.01      
    'latticeAdjustmentDelaunayToleranceF', qParams.tolerance2, ...     % This value is too high and is chosen to reduce the compute time. For production work, either do not set, or set to equal or lower than 0.001
    'marginF', qParams.marginF ...                                     % How much larger lattice to generate so as to minimize artifacts in cone spacing near the edges. For production work, do not set this option.
);
toc

% Save the mosaic for later analysis
if (qParams.saveMosaic)
    save(qParams.mosaicFileName, 'theHexMosaic', '-v7.3');
end

%% Print some grid info
theHexMosaic.displayInfo();

%% Visualize the mosaic, showing both the light collecting area (inner segment) and the geometric area
visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', visualizedAperture, ...
    'apertureShape', 'disks', ...
    'panelPosition', [1 1], 'generateNewFigure', true);

%% Visualize the mosaic with the theoretical cone density plot overlayed
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', visualizedAperture, ...
    'apertureShape', 'disks', ...
    'overlayConeDensityContour', 'theoretical', ...
    'coneDensityContourLevels', (180:20:250)*1000, ...    % cones/mm^2
    'panelPosition', [1 1], 'generateNewFigure', true);

%% Visualize the mosaic with the actual cone density plot overlayed
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', visualizedAperture, ...
    'apertureShape', 'disks', ...
    'overlayConeDensityContour', 'measured', ...
    'coneDensityContourLevels', (180:20:250)*1000, ...      % cones/mm^2
    'panelPosition', [2 1], 'generateNewFigure', false);

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