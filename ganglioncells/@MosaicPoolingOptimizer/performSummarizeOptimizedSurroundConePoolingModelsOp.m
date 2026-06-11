function performSummarizeOptimizedSurroundConePoolingModelsOp(mosaicEccsForSummaryStatistics, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;

    mosaicsData = cell(1, numel(mosaicEccsForSummaryStatistics));

    for iMosaic = 1:numel(mosaicEccsForSummaryStatistics)
        % Get mosaic params
        mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEccsForSummaryStatistics(iMosaic));

        % Generate the mosaic filename
        [mosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
                'mosaicParams', mosaicParams);

        % Load the generated center-only connected mRGCmosaic
        load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

        % Ask the user which optics were used for computing the input cone
        % mosaic STF responses
        opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(...
            mosaicParams, ...
            'opticsChoice', 'n' ...  % Choose native optics
            );

        % Ask the user which H1 cell index was used to optimize the RF
        % surround pooling model models
        [~, gridSamplingScheme,  optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

        % Instantiate a MosaicPoolingOptimizer object with the center-connected
        % mRGC mosaic
        fprintf('\nInstantiating MosaicPoolingOptimizer for mosaic %d/%d....\n', iMosaic, numel(mosaicEccsForSummaryStatistics));
        theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'samplingScheme', gridSamplingScheme, ...
            'generateSamplingGrids', true, ...
            'visualizeSamplingGrids', false);
         fprintf('\nDone instantiating MosaicPoolingOptimizer!\n');

        gridNodesData = cell(1, theMosaicPoolingOptimizerOBJ.gridNodesNum);

        for gridNodeIndex = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum

            LconeRGCindex = theMosaicPoolingOptimizerOBJ.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
            MconeRGCindex = theMosaicPoolingOptimizerOBJ.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);

            optimizedRGCpoolingObjectsFileNameForThisNode = ...
                strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

            fprintf('Loading optimized L- and M-cone RF compute structs from %s\n', optimizedRGCpoolingObjectsFileNameForThisNode);
            
            load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
                'theLconeRFcomputeStruct', ...
                'theMconeRFcomputeStruct');

            s = struct();
            s.LconeRGCposDegs = theMosaicPoolingOptimizerOBJ.theRGCMosaic.rgcRFpositionsDegs(LconeRGCindex,:);
            s.LconeRFcomputeStruct = theLconeRFcomputeStruct;
            s.MconeRGCposDegs = theMosaicPoolingOptimizerOBJ.theRGCMosaic.rgcRFpositionsDegs(MconeRGCindex,:);
            s.MconeRFcomputeStruct = theMconeRFcomputeStruct;

            gridNodesData{gridNodeIndex} = s;
        end % iNode

        mosaicsData{iMosaic} = gridNodesData;
    end % iMosaic

    [hFigSCRcRatioModelFitHistogram, hFigSCIntSensRatioModelFitHistogram] = plotSummarySata(mosaicsData);

    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, 'SCRcRatioModelFitHistogram.pdf');
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigSCRcRatioModelFitHistogram, 300);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, 'SCIntSensRatioModelFitHistogram.pdf');
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFigSCIntSensRatioModelFitHistogram, 300);

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end


function [hFig1, hFig2] = plotSummarySata(mosaicsData)

    achievedSurroundCenterRadiusRatio = nan(numel(mosaicsData), 200);
    achievedSurroundCenterIntSensRatio = achievedSurroundCenterRadiusRatio;
    targetSurroundCenterRadiusRatio = achievedSurroundCenterRadiusRatio;
    targetSurroundCenterIntSensRatio = achievedSurroundCenterRadiusRatio;

    for iMosaic = 1:numel(mosaicsData)
        gridNodesData = mosaicsData{iMosaic};
        for iGridNode = 1:numel(gridNodesData)
            s = gridNodesData{iGridNode};
            targetVisualSTFparams = s.LconeRFcomputeStruct.theTargetSTFparams;
            theFinalSTFdata = s.LconeRFcomputeStruct.theAchievedSTFdata;
            achievedSurroundCenterRadiusRatio(iMosaic, iGridNode) = theFinalSTFdata.fittedDoGModelRsRcRatio;
            achievedSurroundCenterIntSensRatio(iMosaic, iGridNode)  = theFinalSTFdata.fittedDoGModelSCIntSensRatio;
            targetSurroundCenterRadiusRatio(iMosaic, iGridNode) = targetVisualSTFparams.surroundToCenterRcRatio;
            targetSurroundCenterIntSensRatio(iMosaic, iGridNode) = targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio;
        end
    end

    RcRatioRange = [2 10]; 
    IntSensRange = [0.2 0.8];
    

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    dY = (IntSensRange(2)-IntSensRange(1))*ff.axisOffsetFactor;
    IntSensRangeOffset = [IntSensRange(1)+dY IntSensRange(2)];

    dY = (RcRatioRange(2)-RcRatioRange(1))*ff.axisOffsetFactor;
    RcRatioRangeOffset = [RcRatioRange(1)+dY RcRatioRange(2)];

    ax = theAxes{1,1};
    plot(ax, RcRatioRange, RcRatioRange, 'k--', 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    scatter(ax, targetSurroundCenterRadiusRatio(:), achievedSurroundCenterRadiusRatio(:), (ff.markerSize-6)^2, ...
        'LineWidth', ff.lineWidth, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5,  'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1]);

    
    % Range and ticks
    axis(ax, 'equal');
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', RcRatioRangeOffset, 'YLim', RcRatioRangeOffset, 'XTick', 0:2:30, 'YTick', 0:2:30);
    xtickangle(ax, 0);

     % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    xlabel(ax, 'target  r', 'FontAngle', ff.axisFontAngle);
    ylabel(ax, 'achieved  r', 'FontAngle', ff.axisFontAngle);
    

    hFig = figure(2); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    ax = theAxes{1,1};

    histogramEdges = RcRatioRange(1):0.2:RcRatioRange(2);
    baseline = 0.0;
    gain = 100;
    [xBar, yBar, xOutline, yOutline] = generateHistogram(histogramEdges, achievedSurroundCenterRadiusRatio, baseline, gain);

    bar(ax, xBar, yBar, 1, ...
        'FaceColor', [0 0 1], 'EdgeColor',[0 0 1], 'EdgeAlpha', 0.0, ...
        'FaceAlpha', 0.3, 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    plot(ax, xOutline, yOutline, 'b-', 'LineWidth', ff.lineWidth);
    plot(ax, mean(targetSurroundCenterRadiusRatio(:), 'omitnan')*[1 1], [0 100], 'k--', 'LineWidth', ff.lineWidth);

    % Range and ticksdoc
    axis(ax, 'square');
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', RcRatioRangeOffset, 'YLim', 0.3*[ff.axisOffsetFactor 1]*100, 'XTick', 0.:1:30, 'YTick', 0:5:100);
    xtickangle(ax, 0);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    ylabel(ax, 'frequency (%)', 'FontAngle', ff.axisFontAngle);
    xlabel(ax, 'achieved  s', 'FontAngle', ff.axisFontAngle);
    hFig1 = hFig;

    hFig = figure(3); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    ax = theAxes{1,1};
    plot(ax, IntSensRange, IntSensRange, 'k--', 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    scatter(ax, targetSurroundCenterIntSensRatio(:), achievedSurroundCenterIntSensRatio(:), (ff.markerSize-6)^2, ...
        'LineWidth', ff.lineWidth, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.1,   'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1]);


    histogramEdges = IntSensRange(1):0.02:IntSensRange(2);
    baseline = IntSensRange(1);
    gain = 0.5;
    [xBar, yBar, xOutline, yOutline] = generateHistogram(histogramEdges, targetSurroundCenterIntSensRatio, baseline, gain);
    [xBar2, yBar2, xOutline2, yOutline2] = generateHistogram(histogramEdges, achievedSurroundCenterIntSensRatio, baseline, gain);

    % bar(ax, xBar, yBar, 1, ...
    %     'FaceColor', [0 0 1], 'EdgeColor',[0 0 1], 'EdgeAlpha', 0.0, ...
    %     'FaceAlpha', 0.3, 'LineWidth', ff.lineWidth);
    plot(ax, xOutline, yOutline, 'b-', 'LineWidth', ff.lineWidth);

    %bar(ax, yBar2, xBar2, 1, ...
    %    'FaceColor', [0 0 1], 'EdgeColor',[0 0 1], 'EdgeAlpha', 0.0, ...
    %    'FaceAlpha', 0.3, 'LineWidth', ff.lineWidth);
    plot(ax, yOutline2, xOutline2, 'b-', 'LineWidth', ff.lineWidth);

    % Range and ticks
    axis(ax, 'equal');
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', IntSensRangeOffset, 'YLim', IntSensRangeOffset, ...
        'XTick', 0.:0.1:2, 'YTick', 0:0.1:2, ...
        'XTickLabel', sprintf('%0.1f\n', 0.:0.1:2), ...
        'YTickLabel', sprintf('%0.1f\n', 0.:0.1:2));
    xtickangle(ax, 0);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    xlabel(ax, 'target  s', 'FontAngle', ff.axisFontAngle);
    ylabel(ax, 'achieved  s', 'FontAngle', ff.axisFontAngle);
    hFig2 = hFig;
end


function [x, y, xx, yy] = generateHistogram(edges, data, baseline, gain)
    data= data(:);
    idx = find(~isnan(data));
    data = data(idx);
    [counts,edges] = histcounts(data, edges);
    countsPercentage = baseline + counts / sum(counts)*gain;
    dd = (edges(2)-edges(1))/2;
    x = edges(1:end-1)+dd*2;
    y = countsPercentage;

    xx = [];
    yy = [];
    for i = 1:numel(x)-1
        xx(numel(xx)+1) = edges(i)+dd;
        yy(numel(yy)+1) = y(i);
        xx(numel(xx)+1) = edges(i+1)+dd;
        yy(numel(yy)+1) = y(i);
    end

end
