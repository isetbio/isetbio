% Demonstrates how to spatially low-pass a mosaic's responses.
%
% Description:
%    Demonstrates how to spatially low-pass a mosaic's responses. You might
%    want to do this before subsampling, and you might want to subsample to
%    speed up some computation.
%
% See Also:
%    coneMosaic.lowPassMosaicResponse, coneMosaic.demosaicedResponses.

% History:
%    xx/xx/16  NPC  ISETBIO Team, 2016
%    08/08/17  dhb  Remove current computations as they were broken and
%                   not being used. Commenting pass.
%    09/14/18  jnm  Formatting

%% Set up scene
scene = sceneCreate;
scene = sceneSet(scene, 'h fov', 1.0);
oi = oiCreate('human');
oi = oiCompute(oi, scene);

%% Set up mosaic and compute isomerizations
cMosaic = coneMosaic();
cMosaic.setSizeToFOV([sceneGet(scene, 'h fov'), sceneGet(scene, 'v fov')]);
cMosaic.noiseFlag = 'none';
isomerizationsMap = cMosaic.compute(oi);

%% Low-pass the isomerizations map
% The space constants are for L, M and S cones respectively.
lowPassSpaceConstantInMicrons = [5 5 10];
[isomerizationsMapLowPassed, Lmap, Mmap, Smap] = ...
    cMosaic.lowPassMosaicResponse(isomerizationsMap, ...
    lowPassSpaceConstantInMicrons);

%% Visualize results
climRange = [min([min(isomerizationsMap(:)), ...
    min(isomerizationsMapLowPassed(:))]), ...
    max([max(isomerizationsMap(:)), max(isomerizationsMapLowPassed(:))])];
hFig = figure(1);
clf;
set(hFig, 'Position', [10 500 1024 500]);
subplotPosVectors = NicePlot.getSubPlotPosVectors('rowsNum', 2, ...
    'colsNum', 3, 'heightMargin', 0.08, 'widthMargin', 0.05, ...
    'leftMargin', 0.03, 'rightMargin', 0.01, 'bottomMargin', 0.03, ...
    'topMargin', 0.03);

% Original scene
subplot('Position', subplotPosVectors(1, 1).v)
imshow(sceneGet(scene, 'RGB'))
axis 'image';
title('input scene');

% Pre-filtered map of cone isomerizations
subplot('Position', subplotPosVectors(1, 2).v)
imagesc(isomerizationsMap);
set(gca, 'CLim', climRange);
axis 'image';
title(sprintf('original response (max = %2.1f)', ...
    max(isomerizationsMap(:))));

% The map after low pass filtering
subplot('Position', subplotPosVectors(1, 3).v)
imagesc(isomerizationsMapLowPassed);
set(gca, 'CLim', climRange);
axis 'image';
title(sprintf('lowpassed response (max = %2.1f)', ...
    max(isomerizationsMapLowPassed(:))));

% Demosaiced low-pass filtered L plane, with visualization of kernel
subplot('Position', subplotPosVectors(2, 1).v)
imagesc(Lmap);
axis 'image';
title(sprintf('L demosaiced response (lowpassed, tau = %2.1f um)', ...
    lowPassSpaceConstantInMicrons(1)));

% Demosaiced low-pass filtered M plane, with visualization of kernel
subplot('Position', subplotPosVectors(2, 2).v)
imagesc(Mmap);
axis 'image';
title(sprintf('M demosaiced response (lowpassed, tau = %2.1f um)', ...
    lowPassSpaceConstantInMicrons(2)));

% Demosaiced low-pass filtered M plane, with visualization of kernel
subplot('Position', subplotPosVectors(2, 3).v)
imagesc(Smap);
axis 'image';
title(sprintf('S demosaiced response (lowpassed, tau = %2.1f um)', ...
    lowPassSpaceConstantInMicrons(3)));

colormap(gray(1024));
