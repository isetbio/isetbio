function visualizeOptimizationProgress(figNo, figTitle, ...
    targetVisualSTFparams, theCurrentSTFdata, ...
    retinalConePoolingParams,  ...
    pooledConeIndicesAndWeights, ...
    rmseSequence)
        

    ff = MSreadyPlot.figureFormat('2x4');
    hFigProgress = figure(figNo); 

    resetFigure = false;
    if (size(rmseSequence,1) == 1) || (rmseSequence(end,1) < rmseSequence(end-1,1))
        resetFigure = true;
    end

    if (resetFigure)
        clf;
        set(hFigProgress, ...
            'Name', figTitle, ...
            'Position', [10 1000 ff.figureSize(1) ff.figureSize(2)], ...
            'Color', [1 1 1]);
    end

    ax = subplot('Position',  ff.subplotPosVectors(1,1).v);
    MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, retinalConePoolingParams, 'cone pooling');

    ax = subplot('Position',  ff.subplotPosVectors(1,2).v);
    MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, theCurrentSTFdata.fittedDoGModelParams, 'DoG fit to STF');

    if (resetFigure)
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


       % Bottom 
       ff.subplotPosVectors(2,1).v
       largeAxisPosition = [ff.subplotPosVectors(2,1).v(1) ff.subplotPosVectors(2,1).v(2) ...
                            0.43 ff.subplotPosVectors(2,1).v(4)*0.85];
       ax = subplot('Position', largeAxisPosition);
       plot(ax, 1:size(rmseSequence,1), rmseSequence(:,2), 'b.-', 'MarkerSize', 20, 'LineWidth', 1.0); hold(ax, 'on');
       plot(ax, 1:size(rmseSequence,1), rmseSequence(:,3), 'r.-', 'MarkerSize', 20, 'LineWidth', 1.0); hold(ax, 'on');
       xlabel(ax, 'iteration no')
       legend(ax, {'Rs/Rc ratio', 'S/C int. sens. ratio'});
       yTicks = -2:0.2:2;
       set(ax, 'FontSize', ff.fontSize, 'XLim', [0 size(rmseSequence,1)+1], 'YLim', [-1.2 1.2], ...
           'YTick', yTicks, 'YTickLabel', sprintf('%2.1f\n', yTicks));
       grid(ax, 'on');
       box(ax, 'off');

       ax = subplot('Position',  ff.subplotPosVectors(2,3).v);
       axis(ax, 'square');
       title(ax, 'RF center')
       ax = subplot('Position',  ff.subplotPosVectors(2,4).v);
       axis(ax, 'square');
       title(ax, 'RF surround')
    end

   drawnow;
end