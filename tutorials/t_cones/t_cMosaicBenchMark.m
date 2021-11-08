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

% Repeat so as to get an accurate estimate of compute time.  Increase
% this number if you are interested in more accurate timing
repeatsNum = 1;

%% Generate test scenes
%
% Can increase FOV and row/col dimensions to get more dramatic timing
% difference effects, but run time for this tutorial will slow down.
fprintf('Computing scenes ...\n');
sceneFOVDegs = 1.0;
% A sinusoid
params.freq = 10;
params.row = 256;
params.col = 256;
spectrum = 400:2:700;
scene1 = sceneCreate('Harmonic', params, spectrum);
%scene1 = sceneCreate('Harmonic', params);
scene1 = sceneSet(scene1, 'fov', sceneFOVDegs);


scene2 = sceneCreate('distortion grid', 512);
scene2 = sceneSet(scene2, 'fov', sceneFOVDegs);


%% Compute optical images
fprintf('Computing optical images ...\n');
oi1 = oiCreate;
oi1 = oiCompute(scene1, oi1);

oi2 = oiCreate;
oi2 = oiCompute(scene2, oi2);

%% Generate the old-style cone mosaic
fprintf('Computing mosaics ...\n');
regenerateMosaic = true;
resourcesDir = fullfile(strrep(isetRootPath, 'isettools', 'tutorials'), 't_cones', 'resources');
if (regenerateMosaic)
    mosaicFOVDegs = sceneFOVDegs*0.8;
    theOldConeMosaic = coneMosaicHex(7, ...
                'fovDegs', mosaicFOVDegs, ...
                'eccBasedConeDensity', true, ...
                'eccBasedConeQuantalEfficiency', true, ...
                'eccBasedMacularPigment', true, ...
                'maxGridAdjustmentIterations', 20, ...
                'integrationTime', 100/1000, ...
                'noiseFlag', 'none' ...
            );
    save(fullfile(resourcesDir, 'coneMosaicHexBenchMark2.mat'), 'theOldConeMosaic');
else
    load(fullfile(resourcesDir, 'coneMosaicHexBenchMark2.mat'), 'theOldConeMosaic');
end
fprintf(' Done. \n'); 
    
    

%% Compute microns per degree
micronsPerDegFromOpticalImageH = 1e6*oiGet(oi1, 'hspatial resolution')/oiGet(oi1, 'hangular resolution');
micronsPerDegFromOpticalImageW = 1e6*oiGet(oi1, 'wspatial resolution')/oiGet(oi1, 'wangular resolution');
micronsPerDegreeApproximation = 0.5*(micronsPerDegFromOpticalImageH + micronsPerDegFromOpticalImageW);
fprintf('Microns per deg from optical image: %f\n', micronsPerDegreeApproximation);
    
    
%% Generate the new cMosaic using cone params (cone positions and ecc-based 
%% properties) matching the old cone mosaic object. 
theNewConeMosaic = cMosaic(...
            'coneData', theOldConeMosaic.coneData(), ... % coneData from old mosaic
            'micronsPerDegree', micronsPerDegreeApproximation, ...
            'integrationTime', 100/1000, ...
            'noiseFlag', 'none');



%% Visualize the 2 mosaics
fprintf('Visualizing mosaics ...\n');
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
            
%% Setup memory profiler
if (benchtestMemoryAllocation)
    profile('off')
    profile('-memory','on')
end

%% Compute cone excitation responses (old mosaic) to the 2 scenes
fprintf('Computing responses using the OLD mosaic ...\n');

tic
for k = 1:repeatsNum
    oldConeMosaicExcitations1Interleaved = theOldConeMosaic.compute(oi1);
    oldConeMosaicExcitations2Interleaved = theOldConeMosaic.compute(oi2);
end
computeTimeOldConeMosaicInterleaved = toc;

tic
for k = 1:repeatsNum
    oldConeMosaicExcitations1NonInterleaved = theOldConeMosaic.compute(oi1);
end
for k = 1:repeatsNum
    oldConeMosaicExcitations2NonInterleaved = theOldConeMosaic.compute(oi2);
end
computeTimeOldConeMosaicNonInterleaved = toc;


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
fprintf('Computing responses using the NEW mosaic ...\n');

tic
for k = 1:repeatsNum
    newConeMosaicExcitations1Interleaved = theNewConeMosaic.compute(oi1);
    newConeMosaicExcitations2Interleaved = theNewConeMosaic.compute(oi2);
end
computeTimeNewConeMosaicInterleaved = toc;

tic
for k = 1:repeatsNum
    newConeMosaicExcitations1NonInterleaved = theNewConeMosaic.compute(oi1);
end
for k = 1:repeatsNum
    newConeMosaicExcitations2NonInterleaved = theNewConeMosaic.compute(oi2);
end
computeTimeNewConeMosaicNonInterleaved = toc;


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
fprintf('Visualizing results ...\n');
activationRange = [min(newConeMosaicExcitations1Interleaved(:)), max(newConeMosaicExcitations1Interleaved(:))];
ax = subplot(2,2,3);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations1Interleaved, ...
            'signalRange', activationRange, ...
            'mapType', 'modulated disks', ...
            'titleForMap',  ' ');
         
ax = subplot(2,2,4);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations1Interleaved, ...
             'activationRange', activationRange, ...
             'plotTitle',  ' ');

%% Compare responses cone-by-cone
% Concatenate interleaved responses of the old mosaic to the two scenes
lConeIndices = find(theOldConeMosaic.pattern == 2);
mConeIndices = find(theOldConeMosaic.pattern == 3);
sConeIndices = find(theOldConeMosaic.pattern == 4);
oldMosaicLconeResponsesInterleaved = oldConeMosaicExcitations1Interleaved(lConeIndices);
oldMosaicLconeResponsesInterleaved = cat(1, oldMosaicLconeResponsesInterleaved, oldConeMosaicExcitations2Interleaved(lConeIndices));
oldMosaicMconeResponsesInterleaved = oldConeMosaicExcitations1Interleaved(mConeIndices);
oldMosaicMconeResponsesInterleaved = cat(1, oldMosaicMconeResponsesInterleaved, oldConeMosaicExcitations2Interleaved(mConeIndices));
oldMosaicSconeResponsesInterleaved = oldConeMosaicExcitations1Interleaved(sConeIndices);
oldMosaicSconeResponsesInterleaved = cat(1, oldMosaicSconeResponsesInterleaved, oldConeMosaicExcitations2Interleaved(sConeIndices));

% Concatenate interleaved responses of the new mosaic to the two scenes
lConeIndices = theNewConeMosaic.lConeIndices;
mConeIndices = theNewConeMosaic.mConeIndices;
sConeIndices = theNewConeMosaic.sConeIndices;
newMosaicLconeResponsesInterleaved = reshape(newConeMosaicExcitations1Interleaved(lConeIndices), [numel(lConeIndices) 1]);
newMosaicLconeResponsesInterleaved = cat(1, newMosaicLconeResponsesInterleaved, reshape(newConeMosaicExcitations2Interleaved(lConeIndices), [numel(lConeIndices) 1]));
newMosaicMconeResponsesInterleaved = reshape(newConeMosaicExcitations1Interleaved(mConeIndices), [numel(mConeIndices) 1]);
newMosaicMconeResponsesInterleaved = cat(1, newMosaicMconeResponsesInterleaved, reshape(newConeMosaicExcitations2Interleaved(mConeIndices), [numel(mConeIndices) 1]));
newMosaicSconeResponsesInterleaved = reshape(newConeMosaicExcitations1Interleaved(sConeIndices), [numel(sConeIndices) 1]);
newMosaicSconeResponsesInterleaved = cat(1, newMosaicSconeResponsesInterleaved, reshape(newConeMosaicExcitations2Interleaved(sConeIndices), [numel(sConeIndices) 1]));

% Concatenate interleaved responses of the old mosaic to the two scenes
lConeIndices = find(theOldConeMosaic.pattern == 2);
mConeIndices = find(theOldConeMosaic.pattern == 3);
sConeIndices = find(theOldConeMosaic.pattern == 4);
oldMosaicLconeResponsesNonInterleaved = oldConeMosaicExcitations1NonInterleaved(lConeIndices);
oldMosaicLconeResponsesNonInterleaved = cat(1, oldMosaicLconeResponsesNonInterleaved, oldConeMosaicExcitations2NonInterleaved(lConeIndices));
oldMosaicMconeResponsesNonInterleaved = oldConeMosaicExcitations1NonInterleaved(mConeIndices);
oldMosaicMconeResponsesNonInterleaved = cat(1, oldMosaicMconeResponsesNonInterleaved, oldConeMosaicExcitations2NonInterleaved(mConeIndices));
oldMosaicSconeResponsesNonInterleaved = oldConeMosaicExcitations1NonInterleaved(sConeIndices);
oldMosaicSconeResponsesNonInterleaved = cat(1, oldMosaicSconeResponsesNonInterleaved, oldConeMosaicExcitations2NonInterleaved(sConeIndices));

% Concatenate interleaved responses of the new mosaic to the two scenes
lConeIndices = theNewConeMosaic.lConeIndices;
mConeIndices = theNewConeMosaic.mConeIndices;
sConeIndices = theNewConeMosaic.sConeIndices;
newMosaicLconeResponsesNonInterleaved = reshape(newConeMosaicExcitations1NonInterleaved(lConeIndices), [numel(lConeIndices) 1]);
newMosaicLconeResponsesNonInterleaved = cat(1, newMosaicLconeResponsesNonInterleaved, reshape(newConeMosaicExcitations2NonInterleaved(lConeIndices), [numel(lConeIndices) 1]));
newMosaicMconeResponsesNonInterleaved = reshape(newConeMosaicExcitations1NonInterleaved(mConeIndices), [numel(mConeIndices) 1]);
newMosaicMconeResponsesNonInterleaved = cat(1, newMosaicMconeResponsesNonInterleaved, reshape(newConeMosaicExcitations2NonInterleaved(mConeIndices), [numel(mConeIndices) 1]));
newMosaicSconeResponsesNonInterleaved = reshape(newConeMosaicExcitations1NonInterleaved(sConeIndices), [numel(sConeIndices) 1]);
newMosaicSconeResponsesNonInterleaved = cat(1, newMosaicSconeResponsesNonInterleaved, reshape(newConeMosaicExcitations2NonInterleaved(sConeIndices), [numel(sConeIndices) 1]));



%% Plot performance charts
hFig = figure(2); clf;
set(hFig, 'Position', [10 500 800 400]);

% The compute time plot
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
title(sprintf('compute (inteleaved) time speed-up: x %2.1f', computeTimes(1)/computeTimes(2)));

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
ylabel(ax, 'compute time (sec)');
title(sprintf('compute (noninteleaved) time speed-up: x %2.1f', computeTimes(1)/computeTimes(2)));


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


%% Plot correspondence of responses
hFig = figure(3); clf;
set(hFig, 'Position', [10 500 1200 400]);
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,4);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicLconeResponsesInterleaved, newMosaicLconeResponsesInterleaved, 'L-cones (interleaved)');
ax1 = subplot(2,3,2);
ax2 = subplot(2,3,5);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicMconeResponsesInterleaved, newMosaicMconeResponsesInterleaved, 'M-cones (interleaved)')
ax1 = subplot(2,3,3);
ax2 = subplot(2,3,6);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicSconeResponsesInterleaved, newMosaicSconeResponsesInterleaved, 'S-cones (interleaved)')

hFig = figure(4); clf;
set(hFig, 'Position', [10 500 1200 400]);
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,4);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicLconeResponsesNonInterleaved, newMosaicLconeResponsesNonInterleaved, 'L-cones (non-interleaved)');
ax1 = subplot(2,3,2);
ax2 = subplot(2,3,5);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicMconeResponsesNonInterleaved, newMosaicMconeResponsesNonInterleaved, 'M-cones (non-interleaved)')
ax1 = subplot(2,3,3);
ax2 = subplot(2,3,6);
plotCorrespondences(ax1, ax2, activationRange, oldMosaicSconeResponsesNonInterleaved, newMosaicSconeResponsesNonInterleaved, 'S-cones (non-interleaved)')



function plotCorrespondences(ax1, ax2, activationRange, oldMosaicResponse, newMosaicResponse, plotTitle)

    plot(ax1,activationRange, activationRange, 'k-'); hold(ax1, 'on');
    plot(ax1,oldMosaicResponse(:), newMosaicResponse(:), 'r.');
    set(ax1, 'XLim', activationRange, 'YLim', activationRange, 'FontSize', 16);
    set(ax1, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', [0.1 0.3 1 3 10 30 100 300 1000]);
    axis(ax1, 'square');
    grid(ax1, 'on');
    xlabel(ax1, 'excitations (@coneMosaicHex)');
    ylabel(ax1, 'excitations (@cMosaic)');
    title(ax1, plotTitle);

    
    plot(ax2,activationRange, 0*activationRange, 'k-'); hold(ax2, 'on');
    plot(ax2,oldMosaicResponse(:), 100*(newMosaicResponse(:)-oldMosaicResponse(:))./oldMosaicResponse(:), 'r.');
    set(ax2, 'XLim', activationRange, 'FontSize', 16);
    set(ax2, 'XScale', 'log', 'YScale', 'linear', 'XTick', [0.1 0.3 1 3 10 30 100 300 1000], 'YTick', 0.5*(-4:1:4),  'YLim', [-1 1]);
    axis(ax2, 'square');
    grid(ax2, 'on');
    xlabel(ax2, 'excitations (@coneMosaicHex)');
    ylabel(ax2, 'error (%)');
end



