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
examinedMosaicFOVs = [0.25 0.5 1 2];
computeTimesOldConeMosaic = zeros(1, numel(examinedMosaicFOVs));
computeTimesNewConeMosaic = zeros(1, numel(examinedMosaicFOVs));

%% Generate the ring rays stimulus
for iScale = numel(examinedMosaicFOVs):-1:1
    scene = sceneCreate('radial lines', 512);
    scene = sceneSet(scene, 'fov', 1.5*examinedMosaicFOVs(iScale));

    %% Compute the optical image
    oi = oiCreate;
    oi = oiCompute(scene, oi);

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
                'noiseFlag', 'none');
            
    %% Compute cone excitation responses for the 2 mosaics
    repeatsNum = 100;
    tic
    for k = 1:repeatsNum
        oldConeMosaicExcitations = theOldConeMosaic.compute(oi);
    end
    computeTimesOldConeMosaic(iScale) = toc;

    tic
    for k = 1:repeatsNum
        newConeMosaicExcitations = theNewConeMosaic.compute(oi);
    end
    computeTimesNewConeMosaic(iScale) = toc;

end
            
  
        
%% Visualize results
hFig = figure(1); clf;
set(hFig, 'Position', [10 500 1000 400]);
ax = subplot(1,2,1);
plot(ax,examinedMosaicFOVs, computeTimesOldConeMosaic, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7]); hold(ax, 'on');
plot(ax,examinedMosaicFOVs, computeTimesNewConeMosaic, 'ro-','LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log', 'YScale', 'log','FontSize', 16);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'}, 'YTick', [1 3 10 30 100 300]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'compute time (sec)');
legend(ax, {'@coneMosaicHex', '@cMosaic'}, 'Location','northwest');

ax = subplot(1,2,2);
plot(ax,examinedMosaicFOVs, computeTimesOldConeMosaic./computeTimesNewConeMosaic, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); hold(ax, 'on');
set(ax, 'XLim', [0.1 10], 'XScale', 'log',  'FontSize', 16, 'XTick', 1:100);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'});
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'speeding-up factor (x)');
