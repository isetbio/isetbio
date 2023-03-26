function [hFigProgress, ff] = visualizeOptimizationProgress(figNo, figTitle, ...
    targetVisualSTFparams, opticsParams, theCurrentSTFdata, ...
    retinalConePoolingParams,  retinalConePoolingModel, ...
    pooledConeIndicesAndWeights, ...
    rmseSequence)
        
    ff = MSreadyPlot.figureFormat('2x4');
    hFigProgress = figure(figNo); 

    resetFigure = false;
    if (isempty(rmseSequence)) || (size(rmseSequence,1) == 1) || (rmseSequence(end,1) <= min(squeeze(rmseSequence(:,1))))
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

   
    
    if (~isempty(rmseSequence))
        % Bottom left (2,1): The Rs/Rc ratio and the S/C ratio
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
        
        % Bottom left (2,2): The total RMSE
        ax = subplot('Position', ff.subplotPosVectors(2,2).v);
        plot(ax, 1:size(rmseSequence,1), rmseSequence(:,1), 'k.-', 'MarkerSize', 20, 'LineWidth', 1.0);
        xlabel(ax, 'iteration no')
        legend(ax, {'RMSE'}, 'Location', 'SouthWest');
        yTicks = [1e-3 1e-2 1e-1 1 10 100 1000];
        set(ax, 'FontSize', ff.fontSize, 'XLim', [0 size(rmseSequence,1)+1], 'YLim', [1e-3 1000], ...
           'YTick', yTicks, 'YTickLabel', {'1e-3', '1e-2' '0.1', '1', '10', '100', '1e3'}, 'YScale', 'log');
        grid(ax, 'on');
        box(ax, 'off');
    end


    if (resetFigure)
        % These are updated only when the MRSE reaches a new minimum so
        % they depict the best optimization so far during the process

        % top left (1,1): the currently-best H1 surround params
        ax = subplot('Position',  ff.subplotPosVectors(1,1).v);

        if (contains(retinalConePoolingModel, 'H1cellIndex1'))
            conePoolingModelName = 'H1 (#1) surround model params)';
        elseif (contains(retinalConePoolingModel, 'H1cellIndex2'))
            conePoolingModelName = 'H1 (#2) surround model params)';
        elseif (contains(retinalConePoolingModel, 'H1cellIndex3'))
            conePoolingModelName = 'H1 (#3) surround model params)';
        elseif (contains(retinalConePoolingModel, 'H1cellIndex4'))
            conePoolingModelName = 'H1 (#4) surround model params)';
        end
        MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, retinalConePoolingParams, conePoolingModelName);

        % top left (1,2): the currently-best visual STF DoG fit params
        ax = subplot('Position',  ff.subplotPosVectors(1,2).v);
        MosaicPoolingOptimizer.visualizeFittedModelParametersAndRanges(ax, theCurrentSTFdata.fittedDoGModelParams, 'DoG fit to STF');

        % Top right (1,3): The currenty best visual STF
        ax = subplot('Position',  ff.subplotPosVectors(1,3).v);

        if (isfield(theCurrentSTFdata, 'visualizedSpatialFrequencyRange'))
            visualizedSpatialFrequencyRange = theCurrentSTFdata.visualizedSpatialFrequencyRange;
        else
            visualizedSpatialFrequencyRange = [];
        end

        MSreadyPlot.renderSTF(ax, ...
           theCurrentSTFdata.spatialFrequencySupport, ...
           theCurrentSTFdata.visualSTF, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.sfHiRes,...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.compositeSTFHiRes, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.centerSTFHiRes, ...
           theCurrentSTFdata.fittedDoGModelToVisualSTF.surroundSTFHiRes, ...
           sprintf('visual STF\n%s, rank:%d, pupil:%2.1fmm', opticsParams.ZernikeDataBase, opticsParams.examinedSubjectRankOrder, opticsParams.pupilDiameterMM), ...
           {'achieved STF', 'fitted DoG STF', 'fitted center STF', 'fitted surround STF'}, ff, ...
           'noYLabel', true, ...
           'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);

       % Top middle panel: correspondence between achieved and desired DoG ratios at current location                    
       ax = subplot('Position',  ff.subplotPosVectors(2,3).v);
       MSreadyPlot.renderPerformance(ax, ...
                 targetVisualSTFparams.surroundToCenterRcRatio, targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio, ...
                 theCurrentSTFdata.fittedDoGModelRsRcRatio, theCurrentSTFdata.fittedDoGModelSCIntSensRatio, ...
                 ff);
       axis(ax, 'square');

    end

    drawnow;
end