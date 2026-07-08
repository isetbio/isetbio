% Demo how to generate a stereo pair of @cMosaic objects
%
% Description:
%    Builds left- and right-eye cMosaics for the same visual-space
%    position, computes their responses to a letter stimulus, and visualizes
%    the ISETBio convention that cone positions and activations are plotted
%    in visual-space-referred coordinates.

% See Also:
%   t_cMosaicBasic

% History:
%    03/22/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.
%    01/07/24               Still runs in ISETCAM branch
%    07/25/25  NPC  Updated script and expanded the comments

%% Initialize
ieInit;
clear;
close all;

% Whether to generate new mosaics from scratch. Keep the default false for
% routine example smoke tests.
generateNewMosaicLattices = false;
if (generateNewMosaicLattices)
    mosaicFOV = [0.3 0.3];
    mosaicEccInVisualSpace = [-1 0];
    stimOffsetRelativeToMosaicCenter = [0.05 0.05];
    % Generate left mosaic
    cmLeft = cMosaic(...
        'sourceLatticeSizeDegs', 1.0, ...               % Bound the online lattice around the ROI
        'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 5, ...                     % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true, ...          % visualize the convergence
        'useParfor', false ...
        );

    % Generate right mosaic
    cmRight = cMosaic(...
        'sourceLatticeSizeDegs', 1.0, ...               % Bound the online lattice around the ROI
        'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 5, ...                     % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true, ...          % visualize the convergence
        'useParfor', false ...
        );

else
    % Mosaic size
    mosaicFOV = [1.2 0.8];

    % Generate a cone mosaic at a negative horizontal eccentricity in visual space
    mosaicEccInVisualSpace = [-4 0];

    % Where to place the stimulus relative to the mosaic's center (offset in visual degrees)
    stimOffsetRelativeToMosaicCenter = [0.2 0.15];

    % Generate the left eye cone mosaic
    cmLeft = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...                % Use the 58-deg wide lattices
        'whichEye', 'left eye', ...                    
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace, ...
        'useParfor', false);

    % Generate the right eye cone mosaic
    cmRight = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...                % Use the 58-deg wide lattices
        'whichEye', 'right eye', ...                      
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace, ...
        'useParfor', false);
end

hFig = figure(1000);
originalFigureUnits = hFig.Units;
hFig.Units = 'pixels';
set(hFig, 'Position', [10 10 1300 1200]);
hFig.Units = originalFigureUnits;

% Descriptive labels for the two cone mosaics
if (mosaicEccInVisualSpace(1) < 0)
    leftMosaicTitle = 'left cMosaic activation (located in nasal retina)';
    rightMosaicTitle = 'right cMosaic activation (located in temporal retina)';
elseif (mosaicEccInVisualSpace(1) > 0)
    leftMosaicTitle = 'left cMosaic activation (located in temporal retina)';
    rightMosaicTitle = 'right cMosaic activation (located in nasal retina)';
else
    leftMosaicTitle = 'left cMosaic activation';
    rightMosaicTitle = 'right cMosaic activation';
end


% Generate the letter ('A') scene
font = fontCreate;
font = fontSet(font,'character','A');
font = fontSet(font, 'size', 30);
stim = sceneCreate('letter',font);


% Compute the optical image of the letter scene using default human optics
oi = oiCreate('human');
oi = oiCompute(oi,stim,'pad value','mean');

% Compute the left mosaic's activation
leftMosaicActivation = cmLeft.compute(oi, ...
    'opticalImagePositionDegs', mosaicEccInVisualSpace + stimOffsetRelativeToMosaicCenter);

% Compute the right mosaic's activation
rightMosaicActivation = cmRight.compute(oi, ...
    'opticalImagePositionDegs', mosaicEccInVisualSpace + stimOffsetRelativeToMosaicCenter);
 
 
%% Visualize the cone mosaics (in visual space referred coordinates) and the ISETBio meridian conventions
hFig = figure(1001);
originalFigureUnits = hFig.Units;
hFig.Units = 'pixels';
set(hFig, 'Position', [10 10 1600 1150]);
hFig.Units = originalFigureUnits;

% Visualize the meridian labeling convention
ax = subplot('Position', [0.25 0.45 0.55 0.6]);
A = imread('MeridianLabelingConvention.png');
image(ax, A);
set(ax, 'XTick', [], 'YTick', []);
axis(ax, 'image');

% Visualize the left mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.1 0.05 0.4 0.35]);
cmLeft.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture', 'geometricArea', ...
    'densityContourOverlay', true, ...
    'densityContourLevels', 1e3*[50 80 100 130 180 250], ...
    'densityContourLevelLabelsDisplay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmLeft.whichEye);

% Visualize the right mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.55 0.05 0.4 0.35]);
cmRight.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture', 'geometricArea', ...
    'densityContourOverlay', true, ...
    'densityContourLevels', 1e3*[50 80 100 130 180 250], ...
    'densityContourLevelLabelsDisplay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmRight.whichEye);


%% Visualize the stimulus in visual space and the cone mosaic activations (also in visual space referred coordinates)
hFig = figure(1002); clf;
originalFigureUnits = hFig.Units;
hFig.Units = 'pixels';
set(hFig, 'Position', [10 10 1850 1090]);
hFig.Units = originalFigureUnits;

% Visualize the optical image of the stimulus in visual space coordinates
% at the top
ax = subplot('Position', [0.3 0.57 0.45 0.4]);
oiSpatialSupportMeters = oiGet(oi, 'spatial support');
oiSpatialSupportXMicrons = squeeze(oiSpatialSupportMeters(1,1:end,1)) * 1e6;
oiSpatialSupportYMicrons = flipud(squeeze(oiSpatialSupportMeters(1:end,1,2))) * 1e6;
oiSpatialSupportXMicrons = oiSpatialSupportXMicrons + (mosaicEccInVisualSpace(1) + stimOffsetRelativeToMosaicCenter(1))*cmLeft.micronsPerDegree;
oiSpatialSupportYMicrons = oiSpatialSupportYMicrons + (mosaicEccInVisualSpace(2) + stimOffsetRelativeToMosaicCenter(2))*cmLeft.micronsPerDegree;
image(ax, oiSpatialSupportXMicrons,oiSpatialSupportYMicrons, oiGet(oi, 'rgb'));
hold(ax, 'on');

% Superimpose the outline of the cone mosaic
xOutline = cmRight.eccentricityMicrons(1) + 0.5*[-1 -1 1 1 -1]*cmRight.sizeMicrons(1);
yOutline = cmRight.eccentricityMicrons(2) + 0.5*[-1 1 1 -1 -1]*cmRight.sizeMicrons(2);

plot(ax, xOutline, yOutline, 'r-', 'LineWidth', 1.5);
hold(ax, 'off');
xlabel(ax, 'visual space (microns)');
ylabel(ax, 'visual space (microns)');
set(ax, 'XLim', 4000*[-1 1], 'YLim', 2000*[-1 1]);
set(ax, 'FontSize', 16);
axis(ax, 'xy');
title(ax, 'optical image of stimulus (in visual space)');

% Visualize the left mosaic activation
ax = subplot('Position', [0.05 0.05 0.45 0.47]);
cmLeft.visualize('figureHandle', hFig, 'axesHandle', ax, ...
      'activation', leftMosaicActivation, ...
      'visualizedConeAperture', 'geometricArea', ...
      'domain', 'microns', ...
      'horizontalActivationColorBar', ~true, ...
      'labelRetinalMeridians', true, ...
      'plotTitle', leftMosaicTitle);
  
% Visualize the right mosaic activation
ax = subplot('Position', [0.53 0.05 0.45 0.47]);
cmRight.visualize('figureHandle', hFig, 'axesHandle', ax, ...
      'activation', rightMosaicActivation, ...
      'visualizedConeAperture', 'geometricArea', ...
      'domain', 'microns', ...
      'horizontalActivationColorBar', ~true, ...
      'labelRetinalMeridians', true, ...
      'plotTitle', rightMosaicTitle);
  
