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
% 10/04/17  npc   Reorganized.

%% Initialize
ieInit; clear; close all;

%% Set mosaic parameters
mosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...                      % Sets underlying pixel spacing; controls the rectangular sampling of the hex mosaic grid
    'eccBasedConeDensity', true, ...                % Whether to have an eccentricity based, spatially - varying density
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 0.15*300, ...         % Radius of S-cone free retina, in microns (here set to 0.15 deg).
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

quality.tolerance1 = 0.5;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01 
quality.tolerance2 = 0.05;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001 
quality.marginF = [];                               % How much larger lattice to generate so as to minimize artifacts in cone spacing near the edges. If empty, a dynamic adjustment of margin is done for mosaics < 1.0 degs

%% Set import/export options
saveMosaic = false;                                 % whether to save the mosaic
loadMosaic = false;                                 % whether to load a previously saved mosaic
saveMosaicPDF = false;                              % whether to save a PDF of the mosaic

%% Set FOVs examined
fovExamined = [0.4];                                % mosaic FOV % [0.2 0.4 0.8 1.58];

for pIndex = 1:numel(fovExamined)
    mosaicFOV = fovExamined(pIndex);
    mosaicParams.fovDegs = mosaicFOV;
    mosaicFileName = sprintf('mosaic%2.2f.mat', mosaicFOV);
    
    if (loadMosaic)
        load(mosaicFileName);
    else
        tic
        %% Generate the mosaic.  This takes a little while.
        theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
            'name', mosaicParams.name, ...
            'fovDegs', mosaicParams.fovDegs, ...
            'eccBasedConeDensity', mosaicParams.eccBasedConeDensity, ...
            'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ... 
            'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...                   
            'spatialDensity', mosaicParams.spatialDensity, ...
            'latticeAdjustmentPositionalToleranceF', quality.tolerance1, ...         
            'latticeAdjustmentDelaunayToleranceF', quality.tolerance2, ... 
            'maxGridAdjustmentIterations', 50, ...
            'marginF', quality.marginF ... 
        );

        % Save the mosaic for later analysis
        if (saveMosaic)
            save(mosaicFileName, 'theHexMosaic', '-v7.3');
        end
    end
    
    %% Print mosaic info
    theHexMosaic.displayInfo();

    %% Visualize the mosaic, showing both the light collecting area (inner segment) and the geometric area
    visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
    theHexMosaic.visualizeGrid(...
        'visualizedConeAperture', visualizedAperture, ...
        'apertureShape', 'disks', ...
        'panelPosition', [1 1], 'generateNewFigure', true);

    %% Visualize the mosaic with the theoretical and measured cone density plots overlayed
    hFig = theHexMosaic.visualizeGrid(...
        'visualizedConeAperture', visualizedAperture, ...
        'apertureShape', 'disks', ...
        'labelConeTypes', false, ...
        'overlayHexMesh', true, ...
        'overlayConeDensityContour', 'theoretical_and_measured', ...
        'coneDensityContourLevels', 1000*[170 190 210 220 230], ...    % cones/mm^2
        'panelPosition', [2 1], 'generateNewFigure', false);

    %% Export PDF
    if (saveMosaicPDF)
        NicePlot.exportFigToPDF(sprintf('%s.pdf',mosaicFileName), hFig, 300);
    end
end % pIndex
    
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