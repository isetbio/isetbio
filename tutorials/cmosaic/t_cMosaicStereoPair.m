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
%    Still runs in ISETCAM branch (7/23).
%

%% Initialize
ieInit;
clear;
close all;

generateNewMosaics = false;
if (generateNewMosaics)
    mosaicFOV = [1 1];
    mosaicEccInVisualSpace = [0 0];
    % Generate left mosaic
    cmLeft = cMosaic(...
        'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...           % visualize the convergence
        );

    % Geounerate right mosaic
    cmRight = cMosaic(...
        'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
        'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
        'eccentricityDegs', mosaicEccInVisualSpace, ... % ECC:  x=0.0 degs, y= 0.0 degs
        'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
        'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
        'visualizeMeshConvergence', ~true ...            % visualize the convergence
        );

else
    mosaicEccInVisualSpace = [-12 0];
    mosaicFOV = [5 3];

    %mosaicEccInVisualSpace = [0 0];

    cmLeft = cMosaic(...
        'whichEye', 'left eye', ...                    
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace);
    cmRight = cMosaic(...
        'whichEye', 'right eye', ...                      
        'sizeDegs', mosaicFOV, ...
        'eccentricityDegs', mosaicEccInVisualSpace);
end


if (mosaicEccInVisualSpace(1) < 0)
    leftMosaicTitle = 'left cMosaic activation (nasal retina)';
    rightMosaicTitle = 'right cMosaic activation (temporal retina)';
elseif (mosaicEccInVisualSpace(1) > 0)
    leftMosaicTitle = 'left cMosaic activation (temporal retina)';
    rightMosaicTitle = 'right cMosaic activation (retina retina)';
else
    leftMosaicTitle = 'left cMosaic activation';
    rightMosaicTitle = 'right cMosaic activation';
end

if (1==2)
hFig = figure(1000);
set(hFig, 'Position', [10 10 1600 1150]);

% Visualize the left mosaic with the corresponding contour density plot
% superimposed.
ax = subplot('Position', [0.05 0.05 0.9 0.35]);
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
ax = subplot('Position', [0.05 0.55 0.9 0.35]);
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

end


% Generate scene
generateLetterSceneStimulus = true;
if (generateLetterSceneStimulus)
    
    fontSize = 30;
    if (mosaicEccInVisualSpace(1) == 0)
        fontSize = 30;
    end
    font = fontCreate;
    font = fontSet(font,'character','A');
    font = fontSet(font, 'size', fontSize);
    stim = sceneCreate('letter',font);
else
    params.freq = 10;
    params.row = 512;
    params.col = 512;
    spectrum = 400:20:700;
    stim = sceneCreate('Harmonic', params, spectrum);

    % Put grating in the lower bottom part of the field
    spatialSupport = sceneGet(stim, 'spatial support');
    spatialSupportX = squeeze(spatialSupport(1,1:end,1));
    spatialSupportY =-squeeze(spatialSupport(1:end,1,2));
    photonsFullImage = sceneGet(stim, 'photons');
    photons = bsxfun(@plus, photonsFullImage * 0, mean(mean(photonsFullImage,1),2));
    idx = find((spatialSupportX >= max(spatialSupportX)*0.3) & (spatialSupportX <= max(spatialSupportX)*0.8));
    idy = find((spatialSupportY >= max(spatialSupportX)*0.2) & (spatialSupportY <= max(spatialSupportY)*0.6));
    photons(idy, idx,:) = photonsFullImage(idy, idx,:);
    stim = sceneSet(stim, 'photons', photons);

end


%sceneFOVDegs = 8;
%stim = sceneSet(stim, 'fov', sceneFOVDegs);


% Compute the optical image
oi = oiCreate('human');
oi = oiCompute(oi,stim,'pad value','mean');

% Compute mosaic activations
leftMosaicActivation = cmLeft.compute(oi, ...
    'opticalImagePositionDegs', mosaicEccInVisualSpace);
rightMosaicActivation = cmRight.compute(oi, ...
    'opticalImagePositionDegs', mosaicEccInVisualSpace);
 
 
%% Visualize mosaics
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


hFig = figure(1002); clf;
set(hFig, 'Position', [10 10 1850 1090]);

% Visualize optical image of the stimulus in visual space coordinates
ax = subplot('Position', [0.3 0.57 0.45 0.4]);
oiSpatialSupportMeters = oiGet(oi, 'spatial support');
oiSpatialSupportXMicrons = squeeze(oiSpatialSupportMeters(1,1:end,1)) * 1e6;
oiSpatialSupportYMicrons = squeeze(oiSpatialSupportMeters(1:end,1,2)) * 1e6;
image(ax, ...
    oiSpatialSupportXMicrons + mosaicEccInVisualSpace(1)*cmLeft.micronsPerDegree, ...
    oiSpatialSupportYMicrons + mosaicEccInVisualSpace(2)*cmLeft.micronsPerDegree, flipud(oiGet(oi, 'rgb')));
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
  
      
         
