% Demo usage of arbitrary oiPosition in new @cMosaic object
%
% Description:
%    Shows basic usage of how to position stimuli at arbitary retinal
%    positions and compute using the cMosaic object
%
% See Also:
%   t_cMosaicOffAxis

% History:
%    04/03/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

% Create scene
font = fontCreate('A', 'Georgia', 20, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',1.5);

% Generate optical image
oi = oiCreate('human'); 
oi = oiCompute(oi,scene);

% Generate off-axis mosaic
center = [10 10];
cmfov = [2 2];
cm = cMosaic('eccentricityDegs',center,'sizeDegs',cmfov);

% Compute the response, placing the OI at the center of the mosaic
noiseFreeMosaicCentered = cm.compute(oi, 'opticalImagePositionDegs', 'mosaic-centered');

% Compute the response, placing the OI at (x=10.3, y=9.5);
noiseFreeOffCentered1 = cm.compute(oi, 'opticalImagePositionDegs', [10.3 9.5]);

% Compute the response, placing the OI at (x=9.5, y=9.5);
noiseFreeOffCentered2 = cm.compute(oi, 'opticalImagePositionDegs', [9.5 9.5]);

% Plot everything
hFig = figure(100);

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

% The mosaic response when the OI is placed at (x=10.3, y=9.5);
ax = subplot(2,2,3);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', noiseFreeOffCentered1, ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'plotTitle', 'OI: at x=10.3, y=9.5');

% The mosaic response when the OI is placed at (x=9.5, y=9.5);
ax = subplot(2,2,4);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', noiseFreeOffCentered2, ...
    'domain', 'degrees', ...
    'crossHairsOnMosaicCenter', true, ...
    'plotTitle', 'OI: at x=9.5, y=9.5');



