% Benchmark the new @cMosaic object vs the old @coneMosaicHex
%
% Description:
%    Shows how to generate a @cMosaic object that is completely equivalent
%    to a @coneMosaicHex, computes the responses of the @cMosaic and
%    of the original @coneMosaicHex to a scene, and finally contrasts
%    responses and computation times for the 2 mosaics.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxis

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate the ring rays stimulus
scene = sceneCreate('radial lines', 512);
scene = sceneSet(scene, 'fov', 2.0);

%% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

%% Generate the old-style cone mosaic with an FOV of 1 x 1 deg
theOldConeMosaic = coneMosaicHex(7, ...
            'fovDegs', 1.0, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', true, ...
            'eccBasedMacularPigment', true, ...
            'maxGridAdjustmentIterations', 50, ...
            'noiseFlag', 'none' ...
        );
    
%% Generate the new cMosaic using cone params from the old cone mosaic object
theNewConeMosaic = cMosaic(...
            'coneData', theOldConeMosaic.coneData(), ...
            'noiseFlag', 'none');

%% Visualize the 2 mosaics
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1100 900]);
ax = subplot(2,2,1); cla(ax);
theOldConeMosaic.visualizeGrid(...
    'axesHandle', ax, ...
    'ticksInMicrons', true, ...
    'visualizedConeAperture',  'geometricArea', ...
    'plotTitle', '@coneMosaicHex');

ax = subplot(2,2,2); cla(ax);
theNewConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'domainVisualizationTicks', struct('x', -200:50:200, 'y', -200:50:200), ...
    'plotTitle', '@cMosaic');
            
%% Compute cone excitation responses for the 2 mosaics
repeatsNum = 100;
tic
for k = 1:repeatsNum
    oldConeMosaicExcitations = theOldConeMosaic.compute(oi);
end
computeTimeOldConeMosaic = toc;

tic
for k = 1:repeatsNum
    newConeMosaicExcitations = theNewConeMosaic.compute(oi);
end
computeTimeNewConeMosaic = toc;
  
        
%% Visualize responses
activationRange = [min(newConeMosaicExcitations(:)), max(newConeMosaicExcitations(:))];
ax = subplot(2,2,3);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations, ...
            'signalRange', activationRange, ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(2,2,4);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations, ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');

%% Compare responses cone-by-cone
% Extract old mosaic responses
lConeIndices = find(theOldConeMosaic.pattern == 2);
mConeIndices = find(theOldConeMosaic.pattern == 3);
sConeIndices = find(theOldConeMosaic.pattern == 4);
oldMosaicLconeResponses = oldConeMosaicExcitations(lConeIndices);
oldMosaicMconeResponses = oldConeMosaicExcitations(mConeIndices);
oldMosaicSconeResponses = oldConeMosaicExcitations(sConeIndices);

% Extract new mosaic responses
newMosaicLConeResponses = newConeMosaicExcitations(theNewConeMosaic.lConeIndices);
newMosaicMConeResponses = newConeMosaicExcitations(theNewConeMosaic.mConeIndices);
newMosaicSConeResponses = newConeMosaicExcitations(theNewConeMosaic.sConeIndices);  

hFig = figure(2); clf;
set(hFig, 'Position', [10 500 1500 400]);
ax = subplot(1,4,1);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicLconeResponses(:), newMosaicLConeResponses(:), 'r.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
ylabel(ax, 'excitations (@cMosaic)');
title('L-cones');

ax = subplot(1,4,2);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicMconeResponses(:), newMosaicMConeResponses(:), 'g.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('M-cones');

ax = subplot(1,4,3);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicSconeResponses(:), newMosaicSConeResponses(:), 'c.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('S-cones');

ax = subplot(1,4,4);
computeTimes = [computeTimeOldConeMosaic computeTimeNewConeMosaic];
b = bar([1 2], computeTimes);
b.FaceColor = 'flat';
b.CData(1,:) = [.4 .4 .4];
b.CData(2,:) = [1 .2 .4];
set(ax, 'XLim', [0.5 2.5], 'YLim', [0 ceil(max(computeTimes)/5)*5], 'FontSize', 16);
set(ax, 'XTick', [1 2], 'XTickLabel', {'@coneMosaicHex', '@cMosaic'});
set(ax, 'YTick', 0:5:ceil(max(computeTimes)/5)*5);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'compute engine');
ylabel(ax, 'compute time (sec)');
title(sprintf('compute speed-up: x %2.1f', computeTimeOldConeMosaic/computeTimeNewConeMosaic));