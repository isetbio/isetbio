% Benchmark the new @cMosaic object vs the old @coneMosaicHex
%
% Description:
%    Shows how to generate a @cMosaic object that is completely equivalent
%    to a @coneMosaicHex, computes the noise-free responses of the @cMosaic and
%    of the original @coneMosaicHex to two scenes, and contrasts the
%    responses, the computation times, and memory allocation between the 2 mosaics.
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

% Set to true to examine memory allocation between old and new mosaics
% If set to true, the compute time measurements may not be as accurate as
% when set to false, because the memory profiler takes resources
benchtestMemoryAllocation = true;

% Repeat 100 times so as to get an accurate estimate of compute time
repeatsNum = 10;

%% Generate test scenes
fontSize = 20;
dotsPerInch = 96;
display = 'LCD-Apple';
sceneFOVDegs = 2.0;
font = fontCreate('A', 'Georgia', fontSize, dotsPerInch);
scene1  = sceneCreate('letter', font, display);
scene1  = sceneSet(scene1,'fov',sceneFOVDegs);

%scene1 = sceneCreate('radial lines', 512);
%scene1 = sceneSet(scene1, 'fov', sceneFOVDegs);

scene2 = sceneCreate('distortion grid', 512);
scene2 = sceneSet(scene2, 'fov', sceneFOVDegs);

%% Compute optical images
oi1 = oiCreate;
oi1 = oiCompute(scene1, oi1);

oi2 = oiCreate;
oi2 = oiCompute(scene2, oi2);

%% Generate the old-style cone mosaic with an FOV of 1 x 1 deg
fprintf('\nGenerating coneMosaicHex ...');
regenerateMosaic = ~true;
resourcesDir = fullfile(strrep(isetRootPath, 'isettools', 'tutorials'), 't_cones', 'resources');
if (regenerateMosaic)
    mosaicFOVDegs = 1.0;
    theOldConeMosaic = coneMosaicHex(7, ...
                'fovDegs', mosaicFOVDegs, ...
                'eccBasedConeDensity', true, ...
                'eccBasedConeQuantalEfficiency', true, ...
                'eccBasedMacularPigment', true, ...
                'maxGridAdjustmentIterations', 20, ...
                'integrationTime', 100/1000, ...
                'noiseFlag', 'none' ...
            );
    save(fullfile(resourcesDir, 'coneMosaicHexBenchMark.mat'), 'theOldConeMosaic');
else
    load(fullfile(resourcesDir, 'coneMosaicHexBenchMark.mat'), 'theOldConeMosaic');
end
fprintf(' Done. \n');

%% Generate the new cMosaic using cone params (cone positions and ecc-based 
%% properties) matching the old cone mosaic object. 
theNewConeMosaic = cMosaic(...
            'coneData', theOldConeMosaic.coneData(), ... % coneData from old mosaic
            'integrationTime', 100/1000, ...
            'noiseFlag', 'none');

fprintf('Visualizing mosaics ... \n');

%% Visualize the 2 mosaics
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1100 900]);
ax = subplot(3,2,1); cla(ax);
theOldConeMosaic.visualizeGrid(...
    'axesHandle', ax, ...
    'ticksInMicrons', true, ...
    'visualizedConeAperture',  'geometricArea', ...
    'plotTitle', '@coneMosaicHex');

ax = subplot(3,2,2); cla(ax);
theNewConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture',  'geometricArea', ...
    'domainVisualizationTicks', struct('x', -200:50:200, 'y', -200:50:200), ...
    'plotTitle', '@cMosaic');
            
%% Setup memory profiler
if (benchtestMemoryAllocation)
    profile('off')
    profile('-memory','on')
end

%% Compute cone excitation responses (old mosaic) to the 2 scenes
fprintf('Computing using @coneMosaicHex\n');
tic
for k = 1:repeatsNum
    oldConeMosaicExcitations1 = theOldConeMosaic.compute(oi1);
    oldConeMosaicExcitations2 = theOldConeMosaic.compute(oi2);
end
computeTimeOldConeMosaic = toc;


%% Find allocated memory
if (benchtestMemoryAllocation)
    p1 = profile('info');
    for fIndex = 1:numel(p1.FunctionTable)
        if (strcmp(p1.FunctionTable(fIndex).FunctionName, 'coneMosaic.compute'))
            memAllocatedOldConeMosaic = p1.FunctionTable(fIndex).TotalMemAllocated;
        end
    end    
    % Restart profiler
    profile('off')
    profile('-memory','on')
end


%% Compute cone excitation responses (new mosaic) to the 2 scenes
fprintf('Computing using @cMosaic\n');
tic
for k = 1:repeatsNum
    newConeMosaicExcitations1 = theNewConeMosaic.compute(oi1);
    newConeMosaicExcitations2 = theNewConeMosaic.compute(oi2);
end
computeTimeNewConeMosaic = toc;


%% Find allocated memory
if (benchtestMemoryAllocation)
    p2 = profile('info');
    for fIndex = 1:numel(p2.FunctionTable)
        if (strcmp(p2.FunctionTable(fIndex).FunctionName, 'cMosaic.compute'))
            memAllocatedNewConeMosaic = p2.FunctionTable(fIndex).TotalMemAllocated;
        end
    end
end


%% Visualize responses
activationRange = [min(newConeMosaicExcitations1(:)), max(newConeMosaicExcitations1(:))];
ax = subplot(3,2,3);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations1, ...
            'signalRange', activationRange, ...
            'visualizedConeAperture', 'geometricArea', ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(3,2,4);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations1, ...
             'visualizedConeAperture',  'geometricArea', ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');

activationRange = [min(newConeMosaicExcitations2(:)), max(newConeMosaicExcitations2(:))];
ax = subplot(3,2,5);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations2, ...
            'signalRange', activationRange, ...
            'visualizedConeAperture', 'geometricArea', ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(3,2,6);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations2, ...
             'visualizedConeAperture',  'geometricArea', ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');
         
%% Compare responses cone-by-cone
% Concatenate responses of the old mosaic to the two scenes
lConeIndices = find(theOldConeMosaic.pattern == 2);
mConeIndices = find(theOldConeMosaic.pattern == 3);
sConeIndices = find(theOldConeMosaic.pattern == 4);
oldMosaicLconeResponses = oldConeMosaicExcitations1(lConeIndices);
oldMosaicLconeResponses = cat(1, oldMosaicLconeResponses, oldConeMosaicExcitations2(lConeIndices));
oldMosaicMconeResponses = oldConeMosaicExcitations1(mConeIndices);
oldMosaicMconeResponses = cat(1, oldMosaicMconeResponses, oldConeMosaicExcitations2(mConeIndices));
oldMosaicSconeResponses = oldConeMosaicExcitations1(sConeIndices);
oldMosaicSconeResponses = cat(1, oldMosaicSconeResponses, oldConeMosaicExcitations2(sConeIndices));

% Concatenate responses of the new mosaic to the two scenes
lConeIndices = theNewConeMosaic.lConeIndices;
mConeIndices = theNewConeMosaic.mConeIndices;
sConeIndices = theNewConeMosaic.sConeIndices;
newMosaicLconeResponses = reshape(newConeMosaicExcitations1(lConeIndices), [numel(lConeIndices) 1]);
newMosaicLconeResponses = cat(1, newMosaicLconeResponses, reshape(newConeMosaicExcitations2(lConeIndices), [numel(lConeIndices) 1]));
newMosaicMconeResponses = reshape(newConeMosaicExcitations1(mConeIndices), [numel(mConeIndices) 1]);
newMosaicMconeResponses = cat(1, newMosaicMconeResponses, reshape(newConeMosaicExcitations2(mConeIndices), [numel(mConeIndices) 1]));
newMosaicSconeResponses = reshape(newConeMosaicExcitations1(sConeIndices), [numel(sConeIndices) 1]);
newMosaicSconeResponses = cat(1, newMosaicSconeResponses, reshape(newConeMosaicExcitations2(sConeIndices), [numel(sConeIndices) 1]));

%% Plot correspondence of responses
hFig = figure(2); clf;
set(hFig, 'Position', [10 500 1200 400]);
% L-cones
ax = subplot(1,3,1);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicLconeResponses(:), newMosaicLconeResponses(:), 'r.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
ylabel(ax, 'excitations (@cMosaic)');
title('L-cones');

% M-cones
ax = subplot(1,3,2);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicMconeResponses(:), newMosaicMconeResponses(:), 'g.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('M-cones');

% S-cones
ax = subplot(1,3,3);
plot(ax,activationRange, activationRange, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicSconeResponses(:), newMosaicSconeResponses(:), 'c.');
set(ax, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'log');
set(ax, 'XTick', [0.1 0.3 1 3 10 30 100 300 1000 3000 10000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('S-cones');


%% Plot performance charts
hFig = figure(3); clf;
set(hFig, 'Position', [10 500 800 400]);

% The compute time plot
ax = subplot(1,2,1);
computeTimes = [computeTimeOldConeMosaic computeTimeNewConeMosaic];
b = bar([1 2], computeTimes);
b.FaceColor = 'flat';
b.CData(1,:) = [.4 .4 .4];
b.CData(2,:) = [1 .2 .4];
set(ax, 'XLim', [0.5 2.5], 'YLim', [0 (ceil(max(computeTimes)/10)+1)*10], 'FontSize', 16);
set(ax, 'XTick', [1 2], 'XTickLabel', {'@coneMosaicHex', '@cMosaic'});
set(ax, 'YTick', 0:10:ceil(max(computeTimes)/10)*10);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'compute engine');
ylabel(ax, 'compute time (sec)');
title(sprintf('compute time speed-up: x %2.1f', computeTimeOldConeMosaic/computeTimeNewConeMosaic));

% The memory allocation plot
if (benchtestMemoryAllocation)
    ax = subplot(1,2,2);
    allocatedGBytes = [memAllocatedOldConeMosaic memAllocatedNewConeMosaic]/(1024*1024*1024);
    b = bar([1 2], allocatedGBytes);
    b.FaceColor = 'flat';
    b.CData(1,:) = [.4 .4 .4];
    b.CData(2,:) = [1 .2 .4];
    set(ax, 'XLim', [0.5 2.5], 'YLim', [0 (ceil(max(allocatedGBytes)/10)+1)*10], 'FontSize', 16);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'@coneMosaicHex', '@cMosaic'});
    set(ax, 'YTick', 0:10:ceil(max(allocatedGBytes)/10)*10);
    axis(ax, 'square');
    grid(ax, 'on');
    xlabel(ax, 'compute engine');
    ylabel(ax, 'allocated memory (GB)');
    title(sprintf('memory reduction: x %2.1f', allocatedGBytes(1)/allocatedGBytes(2)));
end
