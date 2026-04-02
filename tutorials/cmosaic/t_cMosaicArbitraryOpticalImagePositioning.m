% Demo usage of arbitrary oiPosition in new @cMosaic object
%
% Description:
%    Shows how to position stimuli at arbitary retinal positions and
%    compute using the cMosaic object
%
% See Also:
%   t_cMosaicOffAxis
%
% History:
%    04/03/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

%%
ieInit;

%% Create scene
font = fontCreate('A', 'Georgia', 20, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',1.5);

%% Generate optical image
oi = oiCreate('human'); 
oi = oiCompute(oi,scene,'pad value','mean');

%% Generate off-axis mosaic

% This takes a while because it starts up the parallel workers

center = [4 4];
cmfov  = [2 2];
cm = cMosaic('eccentricityDegs',center,'sizeDegs',cmfov);

% Compute the response, placing the OI at the center of the mosaic
noiseFreeMosaicCentered = cm.compute(oi, 'opticalImagePositionDegs', 'mosaic-centered');

% Compute the response, placing the OI at a different location
noiseFreeOffCentered1 = cm.compute(oi, 'opticalImagePositionDegs', [4.2 3.9]);

% Compute the response, placing the OI at yet a different location
noiseFreeOffCentered2 = cm.compute(oi, 'opticalImagePositionDegs', [3.7 4.3]);

%% Plot everything
hFig = ieFigure;

% The mosaic
ax = subplot(2,2,1);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'plotTitle', 'cone mosaic');

% The mosaic response when the OI is mosaic-centered
ax = subplot(2,2,2);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', noiseFreeMosaicCentered, ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'plotTitle', 'OI: mosaic-centered');

% The mosaic response when the OI is displaced 
ax = subplot(2,2,3);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', noiseFreeOffCentered1, ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'plotTitle', 'OI');

% The mosaic response when the OI is displaced
ax = subplot(2,2,4);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', noiseFreeOffCentered2, ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'plotTitle', 'OI');



