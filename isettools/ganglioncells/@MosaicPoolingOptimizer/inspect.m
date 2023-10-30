function inspect(obj, gridNodeIndex, opticsParams, optimizedRGCpoolingObjectsFileName, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;

    % Optimized RGCpooling object filename
    optimizedRGCpoolingObjectsFileNameForThisNode = ...
        strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

    fprintf('Loading optimized L- and M-cone RF compute structs from %s\n', optimizedRGCpoolingObjectsFileNameForThisNode);
    load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
        'theLconeRFcomputeStruct', ...
        'theMconeRFcomputeStruct');

    LconeRGCindex = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
    fprintf(2,'The inspected L-cone center RGC index is: %d\n', LconeRGCindex);
    figNo = 10000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d L-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), LconeRGCindex);
    pdfFilename = strrep(optimizedRGCpoolingObjectsFileNameForThisNode, '.mat', '_Lcone');

    % Generate the visualization cache
    xSupport = [];
    ySupport = []; 
    centerSubregionContourSamples = 32;
    contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';
    obj.theRGCMosaic.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);

    inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
        opticsParams, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(LconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theLconeRFcomputeStruct, LconeRGCindex, ...
        obj.theRGCMosaic.visualizationCache.rfCenterContourData{LconeRGCindex}, ...
        tickSeparationArcMin, ...
        normalizedPeakSurroundSensitivity, visualizedSpatialFrequencyRange, ...
        gridlessLineWeightingFuncions);

    MconeRGCindex = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
    fprintf(2,'The inspected M-cone center RGC index is: %d\n', MconeRGCindex);
    figNo = 20000 + gridNodeIndex;
    figTitle = sprintf('grid no %d of %d M-cone center RGC %d', ...
        gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), MconeRGCindex);
    pdfFilename = strrep(optimizedRGCpoolingObjectsFileNameForThisNode, '.mat', '_Mcone');

    inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
        opticsParams, ...
        obj.theRGCMosaic.rgcRFpositionsDegs(MconeRGCindex,:), ...
        obj.theRGCMosaic.inputConeMosaic, ...
        theMconeRFcomputeStruct, MconeRGCindex, ...
        obj.theRGCMosaic.visualizationCache.rfCenterContourData{MconeRGCindex}, ...
        tickSeparationArcMin, ...
        normalizedPeakSurroundSensitivity, visualizedSpatialFrequencyRange, ...
        gridlessLineWeightingFuncions);

end

function inspectConeSpecificRFcomputeStruct(figNo, figTitle, pdfFilename, ...
    opticsParams, rgcRFposDegs, inputConeMosaic, theConeSpecificRFcomputeStruct, theRGCindex, ...
    theRFcenterContourData, tickSeparationArcMin, normalizedPeakSurroundSensitivity, visualizedSpatialFrequencyRange, ...
    gridlessLineWeightingFuncions)


    % Retrieve the saved data
    targetVisualSTFparams = theConeSpecificRFcomputeStruct.theTargetSTFparams;
    theFinalSTFdata = theConeSpecificRFcomputeStruct.theAchievedSTFdata;
    theFinalSTFdata.visualizedSpatialFrequencyRange = visualizedSpatialFrequencyRange;
    retinalConePoolingParams = theConeSpecificRFcomputeStruct.retinalConePoolingParams;
    retinalConePoolingModel = theConeSpecificRFcomputeStruct.modelConstants.retinalConePoolingModel;
    theFinalPooledConeIndicesAndWeights = theConeSpecificRFcomputeStruct.theFinalPooledConeIndicesAndWeights;
    rmseSequence = [];

    [hFig, ff] = MosaicPoolingOptimizer.visualizeOptimizationProgress(...
                figNo, figTitle, ...
                targetVisualSTFparams, opticsParams, ...
                theFinalSTFdata, retinalConePoolingParams, retinalConePoolingModel, ...
                theFinalPooledConeIndicesAndWeights, ...
                rmseSequence);

    figure(hFig);

    % Replace the DoG params with a correspondence between H1 cell data and
    % our surround params
    ax = subplot('Position', ff.subplotPosVectors(1,2).v);
    MSreadyPlot.renderFittedH1paramsCorrespondenceToPackerDaceyData(ax, ...
        retinalConePoolingParams, ff);
    
    % Add the RF center cone pooling map
    spatialSupportRangeArcMin = tickSeparationArcMin*4;
    
    ax = subplot('Position',  ff.subplotPosVectors(2,1).v);
    centerLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(ax, ...
        inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.centerConeIndices, ...
        theFinalPooledConeIndicesAndWeights.centerConeWeights, ...
        'withFigureFormat', ff, ...
        'resetAxes', true, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', sprintf('RF center (RGC index: %d)', theRGCindex));

    
    % Add the RF surround cone pooling map
    ax = subplot('Position',  ff.subplotPosVectors(2,2).v);

    surroundLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(ax, ...
        inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.surroundConeIndices, ...
        theFinalPooledConeIndicesAndWeights.surroundConeWeights, ...
        'withFigureFormat', ff, ...
        'resetAxes', true, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', 'RF surround', ...
        'noYLabel', true);


    % Add the line weighting functions
    % Visualized sensitivity range
    %sensitivityRange(1) = -1.05*max([max(surroundLineWeightingFunctions.x.amplitude(:)) max(surroundLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(2) =  1.0*max([max(centerLineWeightingFunctions.x.amplitude(:)) max(centerLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(1) = -normalizedPeakSurroundSensitivity*sensitivityRange(2);

    centerLineWeightingFunctions.x.amplitude = centerLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    centerLineWeightingFunctions.y.amplitude = centerLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.x.amplitude = surroundLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.y.amplitude = surroundLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    sensitivityRange = sensitivityRange / max(sensitivityRange);

    ax = subplot('Position',  ff.subplotPosVectors(1,4).v);
    mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.x, surroundLineWeightingFunctions.x, ...
        sensitivityRange, 'x', ...
        'withFigureFormat', ff, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', 'line weighting functions, X', ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'gridless',gridlessLineWeightingFuncions);

    ax = subplot('Position',  ff.subplotPosVectors(2,4).v);
    mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.y, surroundLineWeightingFunctions.y, ...
        sensitivityRange, 'y', ...
        'withFigureFormat', ff, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', 'line weighting functions, Y', ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'gridless',gridlessLineWeightingFuncions);


    pdfFilename = strrep(pdfFilename, 'intermediateFiles', 'pdfs');
    pdfFilename = sprintf('%s.pdf', pdfFilename);
    NicePlot.exportFigToPDF(pdfFilename, hFig, 300);
    close(hFig);


    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';


    % The achieved residual vector (shape metrics : comparison to Croner & Kaplan)
    pdfFileName = 'OptimizedSurround_ResidualsVector';

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

   
    MSreadyPlot.renderPerformance(ax, ...
                 targetVisualSTFparams.surroundToCenterRcRatio, targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio, ...
                 theFinalSTFdata.fittedDoGModelRsRcRatio, theFinalSTFdata.fittedDoGModelSCIntSensRatio, ...
                 ff, ...
                 'noTitle', true);
    axis(ax, 'square');

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('Unknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The visual STF
    pdfFileName = 'OptimizedSurround_VisualSTF';

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

    theLegends = {'achieved STF', 'fitted DoG STF', 'fitted center STF', 'fitted surround STF'};
    theLegends = {};

    % Compute typical CronerKaplan STF
    typicalCronerKaplanSTF = typicalCronerKaplanSTFAdaptedForCurrentRcDegs(theFinalSTFdata, targetVisualSTFparams);

    MSreadyPlot.renderSTF(ax, ...
           theFinalSTFdata.spatialFrequencySupport, ...
           theFinalSTFdata.visualSTF, ...
           theFinalSTFdata.fittedDoGModelToVisualSTF.sfHiRes,...
           theFinalSTFdata.fittedDoGModelToVisualSTF.compositeSTFHiRes, ...
           theFinalSTFdata.fittedDoGModelToVisualSTF.centerSTFHiRes, ...
           theFinalSTFdata.fittedDoGModelToVisualSTF.surroundSTFHiRes, ...
           '', ...
           theLegends, ff, ...
           'noYLabel', false, ...
           'visualizedSpatialFrequencyRange', theFinalSTFdata.visualizedSpatialFrequencyRange, ...
           'visualizeCronerKaplanTypicalSTF', typicalCronerKaplanSTF);

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The correspondence to H1 cells
    pdfFileName = 'OptimizedSurround_CorrespondenceToPackerDaceyH1cells';

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};
    MSreadyPlot.renderFittedH1paramsCorrespondenceToPackerDaceyData(ax, ...
        retinalConePoolingParams, ff, ...
        'noLegends', true, ...
        'resetAxes', false);

    
    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);



    % The RF center cone pooling weights
    pdfFileName = 'OptimizedSurround_RFcenterWeightsMap';

    hFig = figure(2); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

    centerLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(ax, ...
        inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.centerConeIndices, ...
        theFinalPooledConeIndicesAndWeights.centerConeWeights, ...
        'withFigureFormat', ff, ...
        'resetAxes', false, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', '');

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

    % The RF surround cone pooling weights
    pdfFileName = 'OptimizedSurround_RFsurroundWeightsMap';

    hFig = figure(3); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

    surroundLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(ax, ...
        inputConeMosaic, ...
        rgcRFposDegs, ...
        theFinalPooledConeIndicesAndWeights.surroundConeIndices, ...
        theFinalPooledConeIndicesAndWeights.surroundConeWeights, ...
        'overlayedSubregionContour', theRFcenterContourData, ...
        'withFigureFormat', ff, ...
        'resetAxes', false, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', '');
    

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);



    % Visualized sensitivity range
    %sensitivityRange(1) = -1.05*max([max(surroundLineWeightingFunctions.x.amplitude(:)) max(surroundLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(2) =  1.0*max([max(centerLineWeightingFunctions.x.amplitude(:)) max(centerLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(1) = -normalizedPeakSurroundSensitivity*sensitivityRange(2);
    centerLineWeightingFunctions.x.amplitude = centerLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    centerLineWeightingFunctions.y.amplitude = centerLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.x.amplitude = surroundLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.y.amplitude = surroundLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    sensitivityRange = sensitivityRange / max(sensitivityRange);

    % The line weighting funciton (x) 
    pdfFileName = 'OptimizedSurround_LineWeightingFunctionX';

    hFig = figure(4); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};

    mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.x, surroundLineWeightingFunctions.x, ...
        sensitivityRange, 'x', ...
        'withFigureFormat', ff, ...
        'resetAxes', false, ...
        'xAxisTickAngleRotationDegs', 0, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', '', ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'gridless',gridlessLineWeightingFuncions);

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The line weighting funciton (x) 
    pdfFileName = 'OptimizedSurround_LineWeightingFunctionY';

    hFig = figure(5); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = theAxes{1,1};


    mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(ax, ...
        centerLineWeightingFunctions.y, surroundLineWeightingFunctions.y, ...
        sensitivityRange, 'y', ...
        'withFigureFormat', ff, ...
        'resetAxes', false, ...
        'xAxisTickAngleRotationDegs', 0, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'plotTitle', '', ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'gridless',gridlessLineWeightingFuncions);

    if (strfind(figTitle, 'L-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_LconeCenter.pdf', pdfFileName));
    elseif (strfind(figTitle, 'M-cone center'))
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('%s_MconeCenter.pdf', pdfFileName));
    else
        error('UNknown center dominance: ''%''.', figTitle);
    end
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end


function typicalCronerKaplanSTF = typicalCronerKaplanSTFAdaptedForCurrentRcDegs(theFinalSTFdata, targetVisualSTFparams)
    [~,idx] = ismember('Kc', theFinalSTFdata.fittedDoGModelParams.names);
    kC = theFinalSTFdata.fittedDoGModelParams.finalValues(idx);
    params(1) = kC;

    [~,idx] = ismember('kS/kC', theFinalSTFdata.fittedDoGModelParams.names);
    KSKcRatio = theFinalSTFdata.fittedDoGModelParams.finalValues(idx);
    params(2) = KSKcRatio;

    [~,idx] = ismember('RsToRc', theFinalSTFdata.fittedDoGModelParams.names);
    RsRcRatio = theFinalSTFdata.fittedDoGModelParams.finalValues(idx);
    params(3) = RsRcRatio;

    [~,idx] = ismember('RcDegs', theFinalSTFdata.fittedDoGModelParams.names);
    RcDegs = theFinalSTFdata.fittedDoGModelParams.finalValues(idx);
    params(4) = RcDegs;

    % Exchange the fitted params(2) and (3) with their target values
    params(3) = targetVisualSTFparams.surroundToCenterRcRatio;
    params(2) = targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio/(params(3)^2);

    DoGSTF = @(params,spatialFrequencySupportCPD)(...
                    abs(params(1)       * ( pi * params(4)^2             * exp(-(pi*params(4)*spatialFrequencySupportCPD).^2) ) - ...
                    params(1)*params(2) * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*spatialFrequencySupportCPD).^2) )));
      
    typicalCronerKaplanSTF = DoGSTF(params,theFinalSTFdata.spatialFrequencySupport);

end
