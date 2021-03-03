% Examine how computation time scales with mosaic size for the new @cMosaic
%
% Description:
%    Examine compute times for @cMosaic vs @coneMosaicHex, as a function of
%    mosaic field of view. This one takes a while to run.
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

scalingFactors = [0.25 0.5 1 2 4 8];
computeTimesOldConeMosaic = zeros(1, numel(scalingFactors));
computeTimesNewConeMosaic = zeros(1, numel(scalingFactors));

%% Generate the ring rays stimulus
for iScale = 1:numel(scalingFactors)
    sFactor = scalingFactors(iScale);
    scene = sceneCreate('radial lines', 512);
    scene = sceneSet(scene, 'fov', 1.5*sFactor);

    %% Compute the optical image
    oi = oiCreate;
    oi = oiCompute(scene, oi);

    %% Generate the old-style cone mosaic with an FOV of 1 x 1 deg
    theOldConeMosaic = coneMosaicHex(7, ...
            'fovDegs', 1.0*sFactor, ...
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
plot(ax,scalingFactors, computeTimesOldConeMosaic, 'ks-', 'LineWidth', 1.5); hold(ax, 'on');
plot(ax,scalingFactors, computeTimesNewConeMosaic, 'rs-','LineWidth', 1.5); hold(ax, 'on');
set(ax, 'XLim', [0.1 5], 'XScale', 'log', 'YLim', [2.5 max(computeTimesOldConeMosaic)]*1.2, 'YScale', 'log','FontSize', 16);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'}, 'YTick', [1 3 10 30 100 300]);
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'compute time (sec)');
legend(ax, {'@coneMosaicHex', '@cMosaic'}, 'Location','northwest');

ax = subplot(1,2,2);
plot(ax,scalingFactors, computeTimesOldConeMosaic./computeTimesNewConeMosaic, 'ks-', 'LineWidth', 1.5); hold(ax, 'on');
set(ax, 'XLim', [0.1 5], 'XScale', 'log',  'FontSize', 16, 'XTick', 1:100);
set(ax, 'XTick', [0.25 0.5 1 2 4 8], 'XTickLabel', {'0.25', '0.5', '1.0', '2.0', '4.0', '8.0'});
axis(ax, 'square');
grid(ax, 'on');
xlabel(ax, 'mosaic FOV (degs)');
ylabel(ax, 'speed-up factor (x)');
