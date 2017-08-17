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
% 08/08/17  NPC   Fixed and cleaned up for updated @coneMosaicHex class

%% Initialize
ieInit; clear; close all;

%% Set mosaic parameters
mosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 5, ...                      % Sets underlying pixel spacing; controls the accuracy of the hex mosaic grid (9 is pretty good, but slow)
    'fovDegs', 0.35, ...                            % FOV in degrees
    'eccBasedConeDensity', true, ...                % Whether to have an eccentricity based, spatially - varying density
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 45, ...               % Radius of S-cone free retina, in microns
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

%% Generate the mosaic.  This takes a little while.
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'name', mosaicParams.name, ...
    'fovDegs', mosaicParams.fovDegs, ...
    'eccBasedConeDensity', mosaicParams.eccBasedConeDensity, ...
    'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ... 
    'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...                   
    'spatialDensity', mosaicParams.spatialDensity, ...
    'latticeAdjustmentPositionalToleranceF', 0.01*2, ...        % For best (but much slower results) this should either not get passed or get set to equal or lower than 0.01      
    'latticeAdjustmentDelaunayToleranceF', 0.001*2 ...          % For best (but much slower results) this should either not get passed or get set to equal or lower than 0.001 
);

%% Print some grid info
theHexMosaic.displayInfo();

%% Visualize the mosaic, showing both the light collecting area (inner segment) and the geometric area
visualizedAperture = 'both'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
theHexMosaic.visualizeGrid('visualizedConeAperture', visualizedAperture, 'generateNewFigure', true);

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

