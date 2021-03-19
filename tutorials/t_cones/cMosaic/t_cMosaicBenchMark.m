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
benchtestMemoryAllocation = ~true;

% Repeat 100 times so as to get an accurate estimate of compute time
repeatsNum = 10;

%% Generate test scenes
% Letter A
fontSize = 20;
dotsPerInch = 96;
display = 'LCD-Apple';
sceneFOVDegs = 2.0;
font = fontCreate('A', 'Georgia', fontSize, dotsPerInch);
scene0  = sceneCreate('letter', font, display);
scene0  = sceneSet(scene0,'fov',sceneFOVDegs);

% A sinusoid
params.freq = 10;
params.row = 900;
params.col = 900;
spectrum = 400:2:700;
%scene1 = sceneCreate('Harmonic', params, spectrum);
scene1 = sceneCreate('Harmonic', params);
scene1 = sceneSet(scene1, 'fov', sceneFOVDegs);

% A distortion grid
scene2 = sceneCreate('distortion grid', 900);
scene2 = sceneSet(scene2, 'fov', sceneFOVDegs);

%% Compute optical images
oi1 = oiCreate;
oi1 = oiCompute(scene1, oi1);

oi2 = oiCreate;
oi2 = oiCompute(scene2, oi2);

%% Compute microns per degree
micronsPerDegFromOpticalImageH = 1e6*oiGet(oi1, 'hspatial resolution')/oiGet(oi1, 'hangular resolution');
micronsPerDegFromOpticalImageW = 1e6*oiGet(oi1, 'wspatial resolution')/oiGet(oi1, 'wangular resolution');
micronsPerDegreeApproximation = 0.5*(micronsPerDegFromOpticalImageH + micronsPerDegFromOpticalImageW);
fprintf('Microns per deg from optical image: %f\n', micronsPerDegreeApproximation);

    
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
            'micronsPerDegree', micronsPerDegreeApproximation, ...
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
    oldConeMosaicExcitations1a = theOldConeMosaic.compute(oi1);
end
for k = 1:repeatsNum
    oldConeMosaicExcitations2a = theOldConeMosaic.compute(oi2);
end
computeTimeOldConeMosaicNonInterleaved = toc;

tic
for k = 1:repeatsNum
    oldConeMosaicExcitations1b = theOldConeMosaic.compute(oi1);
    oldConeMosaicExcitations2b = theOldConeMosaic.compute(oi2);
end
computeTimeOldConeMosaicInterleaved = toc;


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
    newConeMosaicExcitations1a = theNewConeMosaic.compute(oi1);
end
for k = 1:repeatsNum
    newConeMosaicExcitations2a = theNewConeMosaic.compute(oi2);
end
computeTimeNewConeMosaicNonInterleaved = toc;

tic
for k = 1:repeatsNum
    newConeMosaicExcitations1b = theNewConeMosaic.compute(oi1);
    newConeMosaicExcitations2b = theNewConeMosaic.compute(oi2);
end
computeTimeNewConeMosaicInterleaved = toc;





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
activationRange = [min(newConeMosaicExcitations1a(:)), max(newConeMosaicExcitations1a(:))];
ax = subplot(3,2,3);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations1a, ...
            'signalRange', activationRange, ...
            'visualizedConeAperture', 'geometricArea', ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(3,2,4);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations1a, ...
             'visualizedConeAperture',  'geometricArea', ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');

activationRange = [min(newConeMosaicExcitations2a(:)), max(newConeMosaicExcitations2a(:))];
ax = subplot(3,2,5);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations2a, ...
            'signalRange', activationRange, ...
            'visualizedConeAperture', 'geometricArea', ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(3,2,6);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations2a, ...
             'visualizedConeAperture',  'geometricArea', ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');
         
%% Compare responses cone-by-cone
% Concatenate responses of the old mosaic to the two scenes
lConeIndices = find(theOldConeMosaic.pattern == 2);
mConeIndices = find(theOldConeMosaic.pattern == 3);
sConeIndices = find(theOldConeMosaic.pattern == 4);
oldMosaicLconeResponsesA = oldConeMosaicExcitations1a(lConeIndices);
oldMosaicLconeResponsesA = cat(1, oldMosaicLconeResponsesA, oldConeMosaicExcitations2a(lConeIndices));
oldMosaicMconeResponsesA = oldConeMosaicExcitations1a(mConeIndices);
oldMosaicMconeResponsesA = cat(1, oldMosaicMconeResponsesA, oldConeMosaicExcitations2a(mConeIndices));
oldMosaicSconeResponsesA = oldConeMosaicExcitations1a(sConeIndices);
oldMosaicSconeResponsesA = cat(1, oldMosaicSconeResponsesA, oldConeMosaicExcitations2a(sConeIndices));

oldMosaicLconeResponsesB = oldConeMosaicExcitations1b(lConeIndices);
oldMosaicLconeResponsesB = cat(1, oldMosaicLconeResponsesB, oldConeMosaicExcitations2b(lConeIndices));
oldMosaicMconeResponsesB = oldConeMosaicExcitations1b(mConeIndices);
oldMosaicMconeResponsesB = cat(1, oldMosaicMconeResponsesB, oldConeMosaicExcitations2b(mConeIndices));
oldMosaicSconeResponsesB = oldConeMosaicExcitations1b(sConeIndices);
oldMosaicSconeResponsesB = cat(1, oldMosaicSconeResponsesB, oldConeMosaicExcitations2b(sConeIndices));

% Concatenate responses of the new mosaic to the two scenes
lConeIndices = theNewConeMosaic.lConeIndices;
mConeIndices = theNewConeMosaic.mConeIndices;
sConeIndices = theNewConeMosaic.sConeIndices;
newMosaicLconeResponsesA = reshape(newConeMosaicExcitations1a(lConeIndices), [numel(lConeIndices) 1]);
newMosaicLconeResponsesA = cat(1, newMosaicLconeResponsesA, reshape(newConeMosaicExcitations2a(lConeIndices), [numel(lConeIndices) 1]));
newMosaicMconeResponsesA = reshape(newConeMosaicExcitations1a(mConeIndices), [numel(mConeIndices) 1]);
newMosaicMconeResponsesA = cat(1, newMosaicMconeResponsesA, reshape(newConeMosaicExcitations2a(mConeIndices), [numel(mConeIndices) 1]));
newMosaicSconeResponsesA = reshape(newConeMosaicExcitations1a(sConeIndices), [numel(sConeIndices) 1]);
newMosaicSconeResponsesA = cat(1, newMosaicSconeResponsesA, reshape(newConeMosaicExcitations2a(sConeIndices), [numel(sConeIndices) 1]));

newMosaicLconeResponsesB= reshape(newConeMosaicExcitations1b(lConeIndices), [numel(lConeIndices) 1]);
newMosaicLconeResponsesB = cat(1, newMosaicLconeResponsesB, reshape(newConeMosaicExcitations2b(lConeIndices), [numel(lConeIndices) 1]));
newMosaicMconeResponsesB = reshape(newConeMosaicExcitations1b(mConeIndices), [numel(mConeIndices) 1]);
newMosaicMconeResponsesB = cat(1, newMosaicMconeResponsesB, reshape(newConeMosaicExcitations2b(mConeIndices), [numel(mConeIndices) 1]));
newMosaicSconeResponsesB = reshape(newConeMosaicExcitations1b(sConeIndices), [numel(sConeIndices) 1]);
newMosaicSconeResponsesB = cat(1, newMosaicSconeResponsesB, reshape(newConeMosaicExcitations2b(sConeIndices), [numel(sConeIndices) 1]));

%% Plot correspondence of responses
hFig = figure(2); clf;
set(hFig, 'Position', [10 500 1200 400]);
% L-cones
ax = subplot(2,3,1);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicLconeResponsesA(:), (newMosaicLconeResponsesA(:)-oldMosaicLconeResponsesA(:))./oldMosaicLconeResponsesA(:)*100, 'r.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
ylabel(ax, 'error (%)');
title('L-cones');

% M-cones
ax = subplot(2,3,2);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicMconeResponsesA(:), (newMosaicMconeResponsesA(:)-oldMosaicMconeResponsesA(:))./oldMosaicMconeResponsesA(:)*100, 'g.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('M-cones');

% S-cones
ax = subplot(2,3,3);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicSconeResponsesA(:), (newMosaicSconeResponsesA(:)-oldMosaicSconeResponsesA(:))./oldMosaicSconeResponsesA(:)*100, 'c.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear');
set(ax, 'XTick', [0.1 0.3 1 3 10 30 100 300 1000 3000 10000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('S-cones');


% L-cones
ax = subplot(2,3,4);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicLconeResponsesB(:), (newMosaicLconeResponsesB(:)-oldMosaicLconeResponsesB(:))./oldMosaicLconeResponsesB(:)*100, 'r.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
ylabel(ax, 'error(%)');
title('L-cones');

% M-cones
ax = subplot(2,3,5);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicMconeResponsesB(:), (newMosaicMconeResponsesB(:)-oldMosaicMconeResponsesB(:))./oldMosaicMconeResponsesB(:)*100, 'g.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('M-cones');

% S-cones
ax = subplot(2,3,6);
plot(ax,activationRange, activationRange*0, 'k-'); hold(ax, 'on');
plot(ax,oldMosaicSconeResponsesB(:), (newMosaicSconeResponsesB(:)-oldMosaicSconeResponsesB(:)) ./ oldMosaicSconeResponsesB(:)*100, 'c.');
set(ax, 'XLim', activationRange, 'YLim', [-1 1], 'FontSize', 16);
set(ax, 'XScale', 'log', 'YScale', 'linear');
set(ax, 'XTick', [0.1 0.3 1 3 10 30 100 300 1000 3000 10000]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'excitations (@coneMosaicHex)');
title('S-cones');


%% Plot performance charts
hFig = figure(3); clf;
set(hFig, 'Position', [10 500 800 400]);

% The inteleaved compute time plot
ax = subplot(1,3,1);
computeTimes = [computeTimeOldConeMosaicInterleaved computeTimeNewConeMosaicInterleaved];
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
title(sprintf('compute time (interleaved) speed-up: x %2.1f', computeTimeOldConeMosaicInterleaved/computeTimeNewConeMosaicInterleaved));

% The non-inteleaved compute time plot
ax = subplot(1,3,2);
computeTimes = [computeTimeOldConeMosaicNonInterleaved computeTimeNewConeMosaicNonInterleaved];
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
ylabel(ax, 'compute time(sec)');
title(sprintf('compute time speed-up (non-inteleaved): x %2.1f', computeTimeOldConeMosaicNonInterleaved/computeTimeNewConeMosaicNonInterleaved));

% The memory allocation plot
if (benchtestMemoryAllocation)
    ax = subplot(1,3,3);
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
