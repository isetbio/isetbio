% Demo how to generate a stereo pair of @cMosaic objects
%
% Description:
%    Shows how to generate a stereo pair of @cMosaic objects, and
%    also demonstrates our convention for labelling visual field & retinal
%    meridians, which can be a bit confusing.
%
% See Also:
%   t_cMosaicBasic

% History:
%    03/22/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

generateNewMosaics = false;
if (generateNewMosaics)
    mosaicFOV = [1 1];

    % Generate left mosaic
    cmLeft = cMosaic(...
        'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...           % visualize the convergence
        );

    % Geounerate right mosaic
    cmRight = cMosaic(...
        'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...            % visualize the convergence
        );

else
    mosaicFOV = [5 3];
    cmLeft = cMosaic(...
        'whichEye', 'left eye', ...                    
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', [0 0]);
    cmRight = cMosaic(...
        'whichEye', 'right eye', ...                      
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', [0 0]);
end

hFig = figure(1000);
set(hFig, 'Position', [10 10 1300 1200]);

% Visualize the left mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.05 0.05 0.9 0.4]);
cmLeft.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture', 'geometricArea', ...
    'densityContourOverlay', true, ...
    'densityContourLevels', 1e3*[20 25 30 35 40 50 60 80 100 130 180 250], ...
    'densityContourLevelLabelsDisplay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmLeft.whichEye);

% Visualize the right mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.05 0.55 0.9 0.4]);
cmRight.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture', 'geometricArea', ...
    'densityContourOverlay', true, ...
    'densityContourLevels', 1e3*[20 25 30 35 40 50 60 80 100 130 180 250], ...
    'densityContourLevelLabelsDisplay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmRight.whichEye);




pause


% Generate scene (sinusoid)
sceneFOVDegs = 4;
params.freq = 10;
params.row = 512;
params.col = 512;
spectrum = 400:20:700;
stim = sceneCreate('Harmonic', params, spectrum);
stim = sceneSet(stim, 'fov', sceneFOVDegs);

% Only put a stimulus in the lower bottom part of the field
spatialSupport = sceneGet(stim, 'spatial support');
spatialSupportX = squeeze(spatialSupport(1,1:end,1));
spatialSupportY =-squeeze(spatialSupport(1:end,1,2));
photonsFullImage = sceneGet(stim, 'photons');
photons = bsxfun(@plus, photonsFullImage * 0, mean(mean(photonsFullImage,1),2));
idx = find((spatialSupportX >= max(spatialSupportX)*0.3) & (spatialSupportX <= max(spatialSupportX)*0.8));
idy = find((spatialSupportY >= max(spatialSupportX)*0.2) & (spatialSupportY <= max(spatialSupportY)*0.6));
photons(idy, idx,:) = photonsFullImage(idy, idx,:);
stim = sceneSet(stim, 'photons', photons);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(stim, oi);

% Compute mosaic activations
leftMosaicActivation = cmLeft.compute(oi);
rightMosaicActivation = cmRight.compute(oi);
 
 
%% Visualize mosaics
hFig = figure(1000);
set(hFig, 'Position', [10 10 1300 1200]);

% Visualize the meridian labeling convention
ax = subplot('Position', [0.25 0.45 0.55 0.6]);
A = imread('MeridianLabelingConvention.png');
image(ax, A);
set(ax, 'XTick', [], 'YTick', []);
axis(ax, 'image');

% Visualize the left mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.1 0.05 0.4 0.4]);
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
ax = subplot('Position', [0.55 0.05 0.4 0.4]);
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


hFig = figure(1001); clf;
set(hFig, 'Position', [10 10 1300 1200]);

% Visualize the optical image
ax = subplot('Position', [0.3 0.57 0.45 0.4]);
oiSpatialSupportMeters = oiGet(oi, 'spatial support');
oiSpatialSupportXMicrons = squeeze(oiSpatialSupportMeters(1,1:end,1)) * 1e6;
oiSpatialSupportYMicrons = squeeze(oiSpatialSupportMeters(1:end,1,2)) * 1e6;
image(ax, oiSpatialSupportXMicrons, oiSpatialSupportYMicrons, flipud(oiGet(oi, 'rgb')));
xlabel(ax, 'retinal space (microns)');
ylabel(ax, 'retinal space (microns)');
set(ax, 'FontSize', 16);
axis(ax, 'xy');
title(ax, 'retinal stimulus');

% Visualize the left mosaic activation
ax = subplot('Position', [0.05 0.05 0.45 0.47]);
cmLeft.visualize('figureHandle', hFig, 'axesHandle', ax, ...
      'activation', leftMosaicActivation, ...
      'visualizedConeAperture', 'geometricArea', ...
      'domain', 'microns', ...
      'horizontalActivationColorBar', true, ...
      'labelRetinalMeridians', true, ...
      'plotTitle', 'left retina activation');
  
% Visualize the right mosaic activation
ax = subplot('Position', [0.53 0.05 0.45 0.47]);
cmRight.visualize('figureHandle', hFig, 'axesHandle', ax, ...
      'activation', rightMosaicActivation, ...
      'visualizedConeAperture', 'geometricArea', ...
      'domain', 'microns', ...
      'horizontalActivationColorBar', true, ...
      'labelRetinalMeridians', true, ...
      'plotTitle', 'right retina activation');
  
      
         
