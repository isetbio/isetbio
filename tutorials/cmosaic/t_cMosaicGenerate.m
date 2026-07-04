% Demo different ways of generating a @cMosaic object
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
    'size degs', [1 1], ...        % A small crop keeps the tutorial quick
    'position degs', [20 -15] ...  % Near the edge of the precomputed lattice
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
    'size degs', [0.25 0.25], ...          % Small enough for interactive use
    'position degs', [1.0 0], ...      % ECC:  x=1.0 degs, y= 0.0 degs
    'compute mesh from scratch', true, ...   % generate mesh on-line, will take some time
    'random seed', 1, ...                   % deterministic tutorial output
    'max mesh iterations', 30 ...           % limit the iterative procedure
    );

%% Visualize in a ieNewGraphWin

hFig = ieNewGraphWin;
cm.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', gca, ...
    'domain', 'degrees', ...
    'plotTitle', 'on-line mesh generation');

drawnow;

%% Experimenting using the renderPatchArray code with the coneMosaicRect

% It will generate a new window
cRect = coneMosaicRect;
coneRectRender(cRect);    % Based on renderPatchArray
xlabel('Position (um)');
ylabel('Position (um)');

%%  Specify the axes

hFig = ieNewGraphWin; thisAxes = gca;
coneRectRender(cRect,'axes handle',thisAxes);
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
