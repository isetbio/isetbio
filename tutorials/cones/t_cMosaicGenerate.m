% Demo different ways of generating a @cMosaic object
%
% Warning:  This script takes a while (few minutes) to run
%
% Description:
%    Shows 3 different ways of generating a @cMosaic object. Also shows how
%    to visualize the generated mosaic.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicStereoPair

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Setup plotting
sv = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', 2, ...
       'heightMargin',  0.09, ...
       'widthMargin',    0.09, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.02);
hFig = figure(1);
set(hFig, 'Position', [10 10 1000 1200]);


%%  Method 1. Generate a @cMosaic object by cropping a region from a large (45x45 deg)
%% precomputed lattice. This is the fastest way to generate a @cMosaic at any eccentricity
cm = cMosaic(...
    'sizeDegs', [4 3], ...            % SIZE: x=4.0 degs, y=3.0 degs
    'eccentricityDegs', [20 -15] ...  % ECC:  x=20 deg, y= -15 deg, near the edge of the precomputed 45x45 mosaic
    );

%% Visualize it (spatial support in degrees)
ax = subplot('Position', sv(1,1).v);
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'degrees', ...
    'plotTitle', 'cropped from large mesh (support: deg)');

%% Visualize it (spatial support in microns)
ax = subplot('Position', sv(1,2).v);
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'plotTitle', 'cropped from large mesh (support: um)');

drawnow;

%%  Method 2. Generate a @cMosic object by generating its mesh from scratch. This can be slow, especially
%% if the mosaic eccentricity is off-axis
cm = cMosaic(...
    'sizeDegs', [0.5 0.5], ...            % SIZE: x=0.5 degs, y=0.5 degs
    'eccentricityDegs', [1.0 0], ...      % ECC:  x=1.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...   % generate mesh on-line, will take some time
    'randomSeed', randi(9999999), ...     % set the random seed, so at to generate a different mosaic each time
    'maxMeshIterations', 80 ...           % stop iterative procedure after this many iterations
    );

%% Visualize it
ax = subplot('Position', sv(2,1).v);
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'degrees', ...
    'plotTitle', 'on-line mesh generation');

drawnow;

%% Method 3. Generate a @coneMosaicHex and its equivalent @cMosaic object
% Generate source mosaic, a @coneMosaicHex object
sourceMosaic = coneMosaicHex(7, 'fovDegs', 0.25);
% Generate equivalent @cMosaic object
cm = cMosaic('coneData', sourceMosaic.coneData());

%% Visualize source @coneMosaicHex mosaic
ax = subplot('Position', sv(3,1).v);
sourceMosaic.visualizeGrid(...
    'axesHandle', ax, ...
    'ticksInMicrons', true, ...
    'visualizedConeAperture',  'lightCollectingArea', ...
    'plotTitle', '@coneMosaicHex (source)');

%% Visualize equivalent @cMosaic
ax = subplot('Position', sv(3,2).v);
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture',  'lightCollectingArea', ...
    'plotTitle', 'source-equivalent @cMosaic');

%%