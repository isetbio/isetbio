% Demo how to generate a stereo pair of @cMosaic objects
%
% Description:
%    Shows how to generate a stereo pair of @cMosaic objects, and
%    also demonstrates our convention for labelling visual field & retinal
%    meridians.
%
% See Also:
%   t_cMosaicBasic

% History:
%    03/22/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

% Generate left mosaic
cmLeft = cMosaic(...
    'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
    'sizeDegs', [1 1], ...                          % SIZE: x=1 degs, y=1 degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', true ...            % visualize the convergence
    );

% Generate right mosaic
cmRight= cMosaic(...
    'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
    'sizeDegs', [1 1], ...                          % SIZE: x=1 degs, y=1 degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', 300, ...                   % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', true ...            % visualize the convergence
    );


%% Visualize mosaics
hFig = figure(1000);
set(hFig, 'Position', [10 10 1300 1030]);

% Visualize the meridian labeling convention
ax = subplot('Position', [0.25 0.45 0.55 0.6]);
A = imread('MeridianLabelingConvention.png');
image(ax, A);
set(ax, 'XTick', [], 'YTick', []);
axis(ax, 'image');

% Visualize the left mosaic
ax = subplot('Position', [0.1 0.05 0.4 0.4]);
cmLeft.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'densityContourOverlay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmLeft.whichEye);

% Visualize the right mosaic
ax = subplot('Position', [0.55 0.05 0.4 0.4]);
cmRight.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'densityContourOverlay', true, ...
    'crossHairsOnFovea', true, ...
    'labelRetinalMeridians', true, ...
    'plotTitle', cmRight.whichEye);
