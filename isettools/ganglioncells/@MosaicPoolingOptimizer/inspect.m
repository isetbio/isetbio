function inspect(obj, gridNodeIndex, optimizedRGCpoolingObjectsFileName, varargin)
    % Parse optional input
%     p = inputParser;
%     p.addParameter('visualize', [], @(x)(isempty(x)||(isscalar(x))));
%     p.parse(varargin{:});

    % Optimized RGCpooling object filename
    optimizedRGCpoolingObjectsFileNameForThisNode = ...
        strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

    fprintf('Loading optimized L- and M-cone RF compute structs from %s\n', optimizedRGCpoolingObjectsFileNameForThisNode);
    load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
        'theLconeRFcomputeStruct', ...
        'theMconeRFcomputeStruct');

    LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
    figNo = 10000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d L-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), LconeRGCindex);
    pdfFilename = strrep(optimizedRGCpoolingObjectsFileNameForThisNode, '.mat', '_Lcone');

    inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(LconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theLconeRFcomputeStruct);

    MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
    figNo = 20000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);
    pdfFilename = strrep(optimizedRGCpoolingObjectsFileNameForThisNode, '.mat', '_Mcone');
    inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(MconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theMconeRFcomputeStruct);

end

function inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
    rgcRFposDegs, inputConeMosaic, theConeSpecificRFcomputeStruct)

    % Retrieve the saved data
    targetVisualSTFparams = theConeSpecificRFcomputeStruct.theTargetSTFparams;
    theFinalSTFdata = theConeSpecificRFcomputeStruct.theAchievedSTFdata;
    retinalConePoolingParams = theConeSpecificRFcomputeStruct.retinalConePoolingParams;
    retinalConePoolingModel = theConeSpecificRFcomputeStruct.modelConstants.retinalConePoolingModel;
    theFinalPooledConeIndicesAndWeights = theConeSpecificRFcomputeStruct.theFinalPooledConeIndicesAndWeights;
    rmseSequence = [];

    [hFig, ff] = MosaicPoolingOptimizer.visualizeOptimizationProgress(...
                figNo, figTitle, ...
                targetVisualSTFparams, ...
                theFinalSTFdata, retinalConePoolingParams, retinalConePoolingModel, ...
                theFinalPooledConeIndicesAndWeights, ...
                rmseSequence);

    figure(hFig);

    % Replace the DoG params with a correspondence between H1 cell data and
    % our surround params
    ax = subplot('Position', ff.subplotPosVectors(1,2).v);
    cla(ax, 'reset');

    idx = find(strcmp(retinalConePoolingParams.names,  'VnVwRatio'));
    fittedModel.NWvolumeRatio = retinalConePoolingParams.finalValues(idx);
    
    idx = find(strcmp(retinalConePoolingParams.names,  'RnRwRatio'));
    fittedModel.RnarrowToRwideRatio = retinalConePoolingParams.finalValues(idx);

    plot(ax, MosaicPoolingOptimizer.PackerDacey2002_H1params.NWvolumeRatios, ...
             MosaicPoolingOptimizer.PackerDacey2002_H1params.RnarrowToRwideRatios, 'ko', ...
             'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5);
    hold(ax, 'on')
    scatter(ax, fittedModel.NWvolumeRatio, fittedModel.RnarrowToRwideRatio, 200, 's', ...
             'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

    xlabel(ax, 'narrow-field/wide-field volume ratio')
    ylabel(ax, 'narrow-field/wide-field radius ratio')
    set(ax, 'XLim', [0 1], 'YLim', [0 0.4], 'XTick', 0:0.2:1, 'YTick', 0:0.1:1, 'TickDir', 'both');
    set(ax, 'FontSize', ff.fontSize)
    axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);
    

    % Add the cone pooling maps
    ax = subplot('Position',  ff.subplotPosVectors(2,3).v);
    cla(ax, 'reset');

    centerLineWeightingFunctions = MosaicPoolingOptimizer.renderConePoolingPlot(ax, inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.centerConeIndices, ...
        theFinalPooledConeIndicesAndWeights.centerConeWeights);
    axis(ax, 'square');
    title(ax, 'RF center');
       
    
    ax = subplot('Position',  ff.subplotPosVectors(2,4).v);
    cla(ax, 'reset');

    surroundLineWeightingFunctions = MosaicPoolingOptimizer.renderConePoolingPlot(ax, inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.surroundConeIndices, ...
        theFinalPooledConeIndicesAndWeights.surroundConeWeights);
    axis(ax, 'square');
    title(ax, 'RF surround')

    % Add the line weighting functions
    sensitivityRange(1) = -1.05*max([max(surroundLineWeightingFunctions.x.amplitude(:)) max(surroundLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(2) = 1.05*max([max(centerLineWeightingFunctions.x.amplitude(:)) max(centerLineWeightingFunctions.y.amplitude(:))]);

    ax = subplot('Position',  ff.subplotPosVectors(2,1).v);
    MosaicPoolingOptimizer.renderConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.x, surroundLineWeightingFunctions.x, ...
        sensitivityRange);
    axis(ax, 'square');
    title(ax, 'line weighting functions (x)')

    ax = subplot('Position',  ff.subplotPosVectors(2,2).v);
    cla(ax, 'reset');

    MosaicPoolingOptimizer.renderConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.y, surroundLineWeightingFunctions.y, ...
        sensitivityRange);
    axis(ax, 'square');
    title(ax, 'line weighting functions (y)')

    
    pdfFilename = strrep(sprintf('%s.pdf',pdfFilename), 'MosaicOptimizerResources', 'MosaicOptimizerPDFs');
    NicePlot.exportFigToPDF(pdfFilename, hFig, 300);

    close(hFig);
end

