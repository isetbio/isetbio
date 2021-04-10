% Examine how computation time scales with mosaic size for the new @cMosaic
%
% Description:
%    Examine compute times for @cMosaic vs @coneMosaicHex, as a function of
%    mosaic field of view. This one takes a long time to run and requires at
%    least 32 GB RAM to run.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicBenchMark
%   t_cMosaicOffAxis

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

% Mosaic sizes to examine
examinedMosaicFOVs = [0.5 1 2 4 8 16];
repeatsNum = 100;
    
computeTimesInterleavedOldConeMosaic = zeros(1, numel(examinedMosaicFOVs));
computeTimesInterleavedNewConeMosaic = zeros(1, numel(examinedMosaicFOVs));
computeTimesNonInterleavedOldConeMosaic = zeros(1, numel(examinedMosaicFOVs));
computeTimesNonInterleavedNewConeMosaic = zeros(1, numel(examinedMosaicFOVs));

%% Generate the ring rays stimulus
for iScale = 1:numel(examinedMosaicFOVs)
    mosaicFOV = examinedMosaicFOVs(iScale);
    [oi1,oi2] = computeOIs(mosaicFOV);
    
    micronsPerDegFromOpticalImageH = 1e6*oiGet(oi1, 'hspatial resolution')/oiGet(oi1, 'hangular resolution');
    micronsPerDegFromOpticalImageW = 1e6*oiGet(oi1, 'wspatial resolution')/oiGet(oi1, 'wangular resolution');
    micronsPerDegreeApproximation = 0.5*(micronsPerDegFromOpticalImageH + micronsPerDegFromOpticalImageW);
    fprintf('Microns per deg from optical image: %f\n', micronsPerDegreeApproximation);

    %% Generate the old-style cone mosaic with an FOV of 1 x 1 deg
    theOldConeMosaic = coneMosaicHex(7, ...
            'fovDegs', examinedMosaicFOVs(iScale), ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', true, ...
            'eccBasedMacularPigment', true, ...
            'maxGridAdjustmentIterations', 20, ...
            'noiseFlag', 'none' ...
        );
    
    %% Generate the new cMosaic using cone params from the old cone mosaic object
    theNewConeMosaic = cMosaic(...
                'coneData', theOldConeMosaic.coneData(), ...
                'micronsPerDegree', micronsPerDegreeApproximation, ...
                'noiseFlag', 'none');
            
    %% Compute cone excitation responses for the 2 mosaics

    tic
    for k = 1:repeatsNum
        oldConeMosaicExcitations = theOldConeMosaic.compute(oi1);
        oldConeMosaicExcitations = theOldConeMosaic.compute(oi2);
    end
    computeTimesInterleavedOldConeMosaic(iScale) = toc;

    tic
    for k = 1:repeatsNum
        newConeMosaicExcitations = theNewConeMosaic.compute(oi1);
        newConeMosaicExcitations = theNewConeMosaic.compute(oi2);
    end
    computeTimesInterleavedNewConeMosaic(iScale) = toc;

    tic
    for k = 1:repeatsNum
        oldConeMosaicExcitations = theOldConeMosaic.compute(oi1);
    end
    for k = 1:repeatsNum
        oldConeMosaicExcitations = theOldConeMosaic.compute(oi2);
    end
    computeTimesNonInterleavedOldConeMosaic(iScale) = toc;
    
    tic
    for k = 1:repeatsNum
        newConeMosaicExcitations = theNewConeMosaic.compute(oi1);
    end
    for k = 1:repeatsNum
        newConeMosaicExcitations = theNewConeMosaic.compute(oi2);
    end
    computeTimesNonInterleavedNewConeMosaic(iScale) = toc;
    
end
            
  
        
%% Visualize results
hFig = figure(1); clf;
set(hFig, 'Position', [10 500 1000 1000]);
ax = subplot(2,2,1);
plot(ax,examinedMosaicFOVs, computeTimesInterleavedOldConeMosaic, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7]); hold(ax, 'on');
plot(ax,examinedMosaicFOVs, computeTimesInterleavedNewConeMosaic, 'ro-','LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log', 'YScale', 'log','FontSize', 16);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'}, 'YTick', [1 3 10 30 100 300]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'compute time (sec)');
legend(ax, {'@coneMosaicHex', '@cMosaic'}, 'Location','northwest');
title(ax, 'interleaved computation');

ax = subplot(2,2,2);
plot(ax,examinedMosaicFOVs, computeTimesInterleavedOldConeMosaic./computeTimesInterleavedNewConeMosaic, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log',  'FontSize', 16, 'XTick', 1:100);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'});
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'speeding-up factor (x)');
title(ax, 'interleaved computation');

ax = subplot(2,2,3);
plot(ax,examinedMosaicFOVs, computeTimesNonInterleavedOldConeMosaic, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7]); hold(ax, 'on');
plot(ax,examinedMosaicFOVs, computeTimesNonInterleavedNewConeMosaic, 'ro-','LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log', 'YScale', 'log','FontSize', 16);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'}, 'YTick', [1 3 10 30 100 300]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'compute time (sec)');
legend(ax, {'@coneMosaicHex', '@cMosaic'}, 'Location','northwest');
title(ax, 'non-interleaved computation');

ax = subplot(2,2,4);
plot(ax,examinedMosaicFOVs, computeTimesNonInterleavedOldConeMosaic./computeTimesNonInterleavedNewConeMosaic, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log',  'FontSize', 16, 'XTick', 1:100);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'});
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'speeding-up factor, non-interleaved (x)');
title(ax, 'non-interleaved computation');

function [oi1,oi2] = computeOIs(mosaicFOV)
    sceneFOV = mosaicFOV * 1.2;
    % A sinusoid
    params.freq = 10;
    params.row = 512;
    params.col = 512;
    spectrum = 400:5:700;
    scene1 = sceneCreate('Harmonic', params, spectrum);
    scene1 = sceneSet(scene1, 'fov', sceneFOV);

    scene2 = sceneCreate('distortion grid', params.row);
    scene2 = sceneSet(scene2, 'fov', sceneFOV);

    oi1 = oiCreate;
    oi1 = oiCompute(scene1, oi1);

    oi2 = oiCreate;
    oi2 = oiCompute(scene2, oi2);
end
