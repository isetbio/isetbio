%% t_coneMosaicHex3
% 
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;

rng('default');
rng(219347);

mosaicParams = struct(...
      'resamplingFactor', 5, ...
                  'size', [32 32], ...
        'spatialDensity', [0 1/3 1/3 1/3]...
    );
    
%% Unit test 1: Mosaic visualization using previous routines
% Generate a rect mosaic 
theRectMosaic = coneMosaic( ...
                   'name', 'the rect mosaic', ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity ...
);

% Generate a hex mosaic using the pattern of the Rect mosaic
theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
                   'name', 'the hex mosaic', ...
                   'size', mosaicParams.size, ...
                'pattern', theRectMosaic.pattern, ...
         'spatialDensity', mosaicParams.spatialDensity ...
);
theHexMosaic.displayInfo();

% Use the standard coneMosaic.plot routine to visualize the mosaics
commandwindow
disp('Hit enter to display the rect and the hex mosaic coneMosaic plotting routines : ');
pause
hFig1 = figure(1);
set(hFig1, 'Position', [10 100 400 400], 'Name', 'Rect mosaic');
hFig2 = figure(2);
set(hFig2, 'Position', [500 100 400 400], 'Name', 'Hex mosaic');
theRectMosaic.plot('cone mosaic', 'hf', hFig1);
theHexMosaic.plot('cone mosaic', 'hf', hFig2);



%% Unit test 2: compare isomerizations between rect and hex mosaics

mosaicParams = struct(...
      'resamplingFactor', 5, ...
                  'size', [48 48], ...
        'spatialDensity', [0 0.62 0.31 0.07], ...
             'noiseFlag', false ...
    );
    
theRectMosaic = coneMosaic( ...
                   'name', 'the rect mosaic', ...
                   'size', mosaicParams.size, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);

theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
                   'name', 'the hex mosaic', ...
                   'size', mosaicParams.size, ...
                'pattern', theRectMosaic.pattern, ...
         'spatialDensity', mosaicParams.spatialDensity, ...
              'noiseFlag', mosaicParams.noiseFlag ...
);

commandwindow
disp('Hit enter to compare isomerizations between rect and hex mosaics (Ring rays scene): ');
pause
%% Ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov', 1.0);
vcAddObject(scene); sceneWindow
    
oi = oiCreate;
oi = oiCompute(scene,oi);   
isomerizationsRect = theRectMosaic.compute(oi,'currentFlag',false);
isomerizationsHex = theHexMosaic.compute(oi,'currentFlag',false);
% render an image that only displays the active cones of the hex mosaic
[activationsHexImage, activationsLMSHexImage] = theHexMosaic.computeActivationImage(isomerizationsHex);

minIsomerization = min([min(isomerizationsRect(:)) min(isomerizationsHex(:))]);
maxIsomerization = max([max(isomerizationsRect(:)) max(isomerizationsHex(:))]);

hFig = figure(100); clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101); clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,1)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,2)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,3)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));




%% Vernier stimulus
commandwindow
disp('Hit enter to compare isomerizations between rect and hex mosaics (Vernier scene): ');
pause

scene = sceneCreate('vernier');
scene.distance = 1;
scene = sceneSet(scene,'fov', 1.0);
vcAddObject(scene); sceneWindow

oi = oiCompute(scene,oi);   
isomerizationsRect = theRectMosaic.compute(oi,'currentFlag',false);
isomerizationsHex = theHexMosaic.compute(oi,'currentFlag',false);
% render an image that only displays the active cones of the hex mosaic
[activationsHexImage, activationsLMSHexImage] = theHexMosaic.computeActivationImage(isomerizationsHex);

minIsomerization = min([min(isomerizationsRect(:)) min(isomerizationsHex(:))]);
maxIsomerization = max([max(isomerizationsRect(:)) max(isomerizationsHex(:))]);

hFig = figure(100); clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101); clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,1)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,2)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,3)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));



%% Slanted bar stimulus
commandwindow
disp('Hit enter to compare isomerizations between rect and hex mosaics (Slanted bar scene): ');
pause

scene = sceneCreate('slanted bar');
scene = sceneSet(scene,'fov', 1.0);
vcAddObject(scene); sceneWindow;

oi = oiCompute(scene,oi);   
isomerizationsRect = theRectMosaic.compute(oi,'currentFlag',false);
isomerizationsHex = theHexMosaic.compute(oi,'currentFlag',false);
% render an image that only displays the active cones of the hex mosaic
[activationsHexImage, activationsLMSHexImage] = theHexMosaic.computeActivationImage(isomerizationsHex);

minIsomerization = min([min(isomerizationsRect(:)) min(isomerizationsHex(:))]);
maxIsomerization = max([max(isomerizationsRect(:)) max(isomerizationsHex(:))]);

hFig = figure(100); clf;
set(hFig, 'Position', [10 10 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.14 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, activationsHexImage);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic isomerizations (all cones)', 'FontSize', 16);

subplot('Position', [0.60 0.03 0.29 0.95]);
imagesc(1:theRectMosaic.cols, 1:theRectMosaic.rows, isomerizationsRect);
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('rect mosaic isomerizations (all cones)', 'FontSize', 16);
colormap(gray(1024));

hFig = figure(101); clf;
set(hFig, 'Position', [10 700 1200 440], 'Color', [1 1 1], 'MenuBar', 'None');
subplot('Position', [0.04 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,1)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic L-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.35 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,2)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic M-cone isomerizations', 'FontSize', 16);

subplot('Position', [0.66 0.03 0.29 0.95]);
imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,3)));
set(gca, 'CLim', [minIsomerization maxIsomerization], 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'
title('hex mosaic S-cone isomerizations', 'FontSize', 16);
colormap(gray(1024));