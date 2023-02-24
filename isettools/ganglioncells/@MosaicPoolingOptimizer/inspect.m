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
    inspectConeSpecificRFcomputeStruct(figNo, figTitle, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(LconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theLconeRFcomputeStruct);

    MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
    figNo = 20000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);
    inspectConeSpecificRFcomputeStruct(figNo, figTitle, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(MconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theMconeRFcomputeStruct);

end

function inspectConeSpecificRFcomputeStruct(figNo, figTitle, ...
    rgcRFposDegs, inputConeMosaic, theConeSpecificRFcomputeStruct)

    % Retrieve the saved data
    targetVisualSTFparams = theConeSpecificRFcomputeStruct.theTargetSTFparams;
    theFinalSTFdata = theConeSpecificRFcomputeStruct.theAchievedSTFdata;
    retinalConePoolingParams = theConeSpecificRFcomputeStruct.retinalConePoolingParams;
    theFinalPooledConeIndicesAndWeights = theConeSpecificRFcomputeStruct.theFinalPooledConeIndicesAndWeights;
    rmseSequence = [];

    hFig = MosaicPoolingOptimizer.visualizeOptimizationProgress(...
                figNo, figTitle, ...
                targetVisualSTFparams, ...
                theFinalSTFdata, retinalConePoolingParams, ...
                theFinalPooledConeIndicesAndWeights, ...
                rmseSequence);


    % Add the cone pooling RF maps
    ff = MSreadyPlot.figureFormat('2x4');
    figure(hFig);
    ax = subplot('Position',  ff.subplotPosVectors(2,3).v);
    
    centerLineWeightingFunctions = MosaicPoolingOptimizer.renderConePoolingPlot(ax, inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.centerConeIndices, ...
        theFinalPooledConeIndicesAndWeights.centerConeWeights);
    axis(ax, 'square');
    title(ax, 'RF center');
       
    
    ax = subplot('Position',  ff.subplotPosVectors(2,4).v);
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
    MosaicPoolingOptimizer.renderConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.y, surroundLineWeightingFunctions.y, ...
        sensitivityRange);
    axis(ax, 'square');
    title(ax, 'line weighting functions (y)')

    NicePlot.exportFigToPDF(sprintf('%s.pdf',figTitle), hFig, 300);
    close(hFig);
end

