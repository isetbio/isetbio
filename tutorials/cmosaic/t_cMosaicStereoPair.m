% Demo how to generate a stereo pair of @cMosaic objects
%
% Description:
%    Shows how to generate a stereo pair of @cMosaic objects, and
%    also demonstrates our convention for labelling visual field & retinal
%    meridians, which can be a bit confusing.
%
%    In ISETBio it was decided early on that all coordinates, even those of cone positions, 
%    are visual space -- referred. This can be a bit confusing, so this tutorial attempts to clarify things. 
%    The first figure depicts the conventions used in ISETBio.
%    The second figure shows an example of responses for a pair of a left
%    and a right eye cone mosaics.
%
%    The top plot of the second figure depicts a stimulus (letter A), positioned, in VISUAL SPACE, 
%    at a negative horizontal eccentricity and at a positive vertical eccentricity. 
%    The red outline depicts the outline of a model cone mosaic, projected in visual space. 
%    The code generates two cone mosaics that correspond to the red outline, 
%    with one mosaic located in the left eye and the other mosaic located in the right eye.
%    The stimulus is projected in both eyes, and its optical image is captured by the two cone mosaics.
%
%    The activations of the the left and right eye cone mosaics are depicted at the bottom two plots. 
%    For the left eye, the stimulus lands in the nasal part of the retina, near the optic disk, 
%    where the cone mosaic model has no cones. In the right eye, the stimulus lands in the temporal part 
%    of the retina. Since the cone mosaics and their activations are depicted in visual space - referred coordinates in ISETBio, 
%    there is no up-down reversal of cone activations with respect to the stimulus. 
%    But the up-down reversal can be evidenced if you note the labeling of the retinal meridians 
%    in these activation plots, which are shown in color. The inferior retina label (in yellow) 
%    appears at positive coordinates, and the superior retinal label (in blue) appears at negative coordinates. 
%    
%    Along the same lines, since the cone mosaic activations are depicted in VISUAL space, 
%    the left and right mosaic activations are not depicted as mirror images of each other. 
%    But you can again note, that the temporal side of the left mosaic (located in the nasal retina), 
%    which is shown in green, appears on the right side, 
%    whereas for the right cone mosaic (located in the temporal retina), the temporal side, 
%    which is shown in red, appears on the left side.

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

% Whether to generate new mosaics from scratch (can take a very long time for large mosaicFOVs)
generateNewMosaicLattices = false;
if (generateNewMosaicLattices)
    mosaicFOV = [1 1];
    mosaicEccInVisualSpace = [0 0];
    % Generate left mosaic
    cmLeft = cMosaic(...
        'sourceLatticeSizeDegs', 64, ...                % Use the 64-deg wide lattices
        'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...           % visualize the convergence
        );

    % Generate right mosaic
    cmRight = cMosaic(...
        'sourceLatticeSizeDegs', 64, ...                % Use the 64-deg wide lattices
        'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...            % visualize the convergence
        );

else
    % Mosaic size
    mosaicFOV = [5 3];

    % Generate a cone mosaic at a negativ e horizontal eccentricity in visual space
    mosaicEccInVisualSpace = [-12 0];

    % Where to place the stimulus relative to the mosaic's center (offset in visual degrees)
    stimOffsetRelativeToMosaicCenter = [1 1];

    % Generate the left eye cone mosaic
    cmLeft = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...                % Use the 58-deg wide lattices
        'whichEye', 'left eye', ...                    
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace);

    % Generate the right eye cone mosaic
    cmRight = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...                % Use the 58-deg wide lattices
        'whichEye', 'right eye', ...                      
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace);
end


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
set(hFig, 'Position', [10 10 1600 1150]);

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
set(hFig, 'Position', [10 10 1850 1090]);

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
  
