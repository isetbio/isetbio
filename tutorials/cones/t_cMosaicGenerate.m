% Demo different ways of generating a @cMosaic object
%
% Warning:  This script takes a while (few minutes) to run
%
% Description:
%    Shows 3 different ways of generating a @cMosaic object. Also shows how
%    to visualize the generated mosaic.
%
% coneMosaicHex - not working in isetcam branch
%
% See Also:
%   t_cMosaicBasic.mlx - for the advanced cone mosaic methods
%   t_cMosaicStereoPair

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;

%%  Method 1. Generate a @cMosaic object 

% This method crops a region from a large (45x45 deg) precomputed
% lattice. This is the fastest way to generate a @cMosaic at any
% eccentricity
cm = cMosaic(...
    'size degs', [4 3], ...            % SIZE: x=4.0 degs, y=3.0 degs
    'position degs', [20 -15] ...  % ECC:  x=20 deg, y= -15 deg, near the edge of the precomputed 45x45 mosaic
    );

%% Visualize it (spatial support in degrees)
cm.visualize(...
    'domain','degrees',...
    'plot title','Support: deg', ...    
    'visualized cone aperture theta samples', 12);

%% Visualize the mosaic (spatial support in microns)

cm.visualize(...
    'domain', 'microns', ...
    'plot title', 'Support: um');
drawnow;

%%  Method 2. Generate a @cMosic object from scratch

% Generating from scratch can be slow, especially
% if the mosaic eccentricity is off-axis
cm = cMosaic(...
    'size degs', [0.5 0.5], ...            % SIZE: x=0.5 degs, y=0.5 degs
    'position degs', [1.0 0], ...      % ECC:  x=1.0 degs, y= 0.0 degs
    'compute mesh from scratch', true, ...   % generate mesh on-line, will take some time
    'random seed', randi(9999999), ...     % set the random seed, so at to generate a different mosaic each time
    'max mesh iterations', 80 ...           % stop iterative procedure after this many iterations
    );

%% Visualize in a ieNewGraphWin

hFig = ieNewGraphWin;
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', gca, ...
    'domain', 'degrees', ...
    'plotTitle', 'on-line mesh generation');

drawnow;

%% Experimenting with using renderPatchArray and the coneMosaicRect


cRect = coneMosaicRect;
%{
% Circular apoerture shape
deltaAngle = 45;
iTheta = (0:deltaAngle:360) / 180 * pi;
coneApertureShape.x = cos(iTheta);
coneApertureShape.y = sin(iTheta);

rfPositions = cRect.coneLocs;
lConeIndices = (cRect.pattern == 2);
mConeIndices = (cRect.pattern == 3);
sConeIndices = (cRect.pattern == 4);
kConeIndices = (cRect.pattern == 1);

% Pattern sample size is the aperture of each cone (square).  In this
% case, they are all the same.
diameter = cRect.patternSampleSize(1);

hFig = ieNewGraphWin; axesHandle = gca;

edgeColor = [0.1, 0.1, 0.1];
lineWidth = 1;
faceAlpha = 1;
edgeAlpha = 1;

faceColors = 1; % 1/4*0.7;
rfCoords = rfPositions(lConeIndices,:);
apertureRadii = ones(size(rfCoords,1),1)*diameter/2;
coneRectRender(axesHandle, coneApertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, faceAlpha, edgeAlpha);

rfCoords = rfPositions(mConeIndices,:);
apertureRadii = ones(size(rfCoords,1),1)*diameter/2;
faceColors = 2; % 2/4*0.7;

coneRectRender(axesHandle, coneApertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, faceAlpha, edgeAlpha);

rfCoords = rfPositions(sConeIndices,:);
apertureRadii = ones(size(rfCoords,1),1)*diameter/2;
faceColors = 3; % 3/4*0.7;
coneRectRender(axesHandle, coneApertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, faceAlpha, edgeAlpha);
%}
hFig = ieNewGraphWin;
thisAxes = gca;
coneRectRender(cRect,thisAxes);
coneRectRender(cRect);
xlabel('Position (um)');
ylabel('Position (um)');


%% Method 3. Generate a @coneMosaicHex and its equivalent @cMosaic object

% coneMosaicHex not working yet.

%{
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
%}