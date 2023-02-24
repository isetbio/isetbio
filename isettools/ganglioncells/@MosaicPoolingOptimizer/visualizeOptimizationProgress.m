function hFigProgress = visualizeOptimizationProgress(figNo, figTitle, ...
    targetVisualSTFparams, theCurrentSTFdata, ...
    retinalConePoolingParams,  ...
    pooledConeIndicesAndWeights, ...
    rmseSequence)
        

    ff = MSreadyPlot.figureFormat('2x4');
    hFigProgress = figure(figNo); 

    resetFigure = false;
    if (size(rmseSequence,1) == 1) || (rmseSequence(end,1) <= min(squeeze(rmseSequence(:,1))))
        resetFigure = true;
    end

    if (resetFigure)
        clf;
        xOffset = 10;
        yOffset = 1000;
        if (mod(figNo,2) == 0)
            xOffset = 1400;
            yOffset = 10;
        end
        set(hFigProgress, ...
            'Name', figTitle, ...
            'Position', [xOffset yOffset ff.figureSize(1) ff.figureSize(2)], ...
            'Color', [1 1 1]);
    end

   
    % Bottom 
    ax = subplot('Position', ff.subplotPosVectors(2,1).v);
    plot(ax, 1:size(rmseSequence,1), rmseSequence(:,2), 'b.-', 'MarkerSize', 20, 'LineWidth', 1.0); hold(ax, 'on');
    plot(ax, 1:size(rmseSequence,1), rmseSequence(:,3), 'r.-', 'MarkerSize', 20, 'LineWidth', 1.0); hold(ax, 'on');
    xlabel(ax, 'iteration no')
    legend(ax, {'Rs/Rc ratio', 'S/C int. sens. ratio'}, 'Location', 'SouthWest');
    yTicks = -2:0.2:2;
    set(ax, 'FontSize', ff.fontSize, 'XLim', [0 size(rmseSequence,1)+1], 'YLim', [-1.2 1.2], ...
       'YTick', yTicks, 'YTickLabel', sprintf('%2.1f\n', yTicks));
    grid(ax, 'on');
    box(ax, 'off');
    
    ax = subplot('Position', ff.subplotPosVectors(2,2).v);
    plot(ax, 1:size(rmseSequence,1), rmseSequence(:,1), 'k.-', 'MarkerSize', 20, 'LineWidth', 1.0);
    xlabel(ax, 'iteration no')
    legend(ax, {'RMSE'}, 'Location', 'SouthWest');
    yTicks = [1e-4 1e-3 1e-2 1e-1 1 10 100];
    set(ax, 'FontSize', ff.fontSize, 'XLim', [0 size(rmseSequence,1)+1], 'YLim', [1e-4 100], ...
       'YTick', yTicks, 'YTickLabel', sprintf('%2.1f\n', yTicks), 'YScale', 'log');
    grid(ax, 'on');
    box(ax, 'off');

    if (resetFigure)
        % These are updated only when the MRSE reaches a new minimum so
        % they depict the best optimization so far during the process
        ax = subplot('Position',  ff.subplotPosVectors(1,1).v);
        MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, retinalConePoolingParams, 'cone pooling');

        ax = subplot('Position',  ff.subplotPosVectors(1,2).v);
        MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, theCurrentSTFdata.fittedDoGModelParams, 'DoG fit to STF');


        ax = subplot('Position',  ff.subplotPosVectors(1,3).v);
        MSreadyPlot.renderSTF(ax, ...
           theCurrentSTFdata.spatialFrequencySupport, ...
           theCurrentSTFdata.visualSTF, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.compositeSTF, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.centerSTF, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.surroundSTF, ...
           '', ...
           {'achieved STF', 'fitted DoG STF', 'fitted center STF', 'fitted surround STF'}, ff, ...
           'noYLabel', true);

       % Top middle panel: correspondence between achieved and desired DoG ratios at current location                    
       ax = subplot('Position',  ff.subplotPosVectors(1,4).v);
       MSreadyPlot.renderPerformance(ax, ...
                 targetVisualSTFparams.surroundToCenterRcRatio, targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio, ...
                 theCurrentSTFdata.fittedDoGModelRsRcRatio, theCurrentSTFdata.fittedDoGModelSCIntSensRatio, ...
                 ff);
       axis(ax, 'square');

       ax = subplot('Position',  ff.subplotPosVectors(2,3).v);
       axis(ax, 'square');
       title(ax, 'RF center')
       ax = subplot('Position',  ff.subplotPosVectors(2,4).v);
       axis(ax, 'square');
       title(ax, 'RF surround')

    end

    drawnow;
end