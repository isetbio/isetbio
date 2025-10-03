function [mRGCindicesToVisualizeSTFsAcrossChromaticities] = findMRGCindicesWithDesiredSurroundConeMix(...
        theComputeReadyMRGCmosaic, targetRangeForSurroundConeMix, targetRangeForCenterMix, ...
        performSurroundAnalysisForConesExclusiveToTheSurround, ...
        rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript)

    % Analyze the surround cone mix for all mRGCs in theComputeReadyMRGCmosaic
    [surroundConeMix, theCenterMajorityConeType, centerConeMix] = surroundConeMixForAllCellsInMosaic(...
            theComputeReadyMRGCmosaic, performSurroundAnalysisForConesExclusiveToTheSurround);
    
    % Determine visualized mRGCindices (only those whose surround cone mix is within the target range)
    mRGCindicesToVisualizeSTFsAcrossChromaticities = find(...
        (surroundConeMix >= targetRangeForSurroundConeMix(1)) & ... 
        (surroundConeMix <= targetRangeForSurroundConeMix(2)) & ...
        (centerConeMix >= targetRangeForCenterMix(1)) & ...
        (centerConeMix <= targetRangeForCenterMix(2)));

    % Visualize the surround cone mix histograms for all L-center and
    % M-center cells in this mRGCmosaic
    pdfFileName = sprintf('SurroundConeMix_eccDegs_%2.1f_%2.1f.pdf', theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(1));

     
    hFig = figure(555); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff, 'figPosition', [576 707]);
    plotTitle = '';
    MSreadyPlot.renderSurroundMixHistograms(theAxes{1,1}, ...
            surroundConeMix(find(theCenterMajorityConeType == cMosaic.LCONE_ID)), ...
            surroundConeMix(find(theCenterMajorityConeType == cMosaic.MCONE_ID)), ...
            plotTitle, ff, ...
            'targetRangeForSurroundConeMix', targetRangeForSurroundConeMix);


    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFig, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileName);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    % Find L-center RGC indices with a surround cone mix in the targetRangeForSurroundConeMix
    idxLconeCenter = find(theCenterMajorityConeType(mRGCindicesToVisualizeSTFsAcrossChromaticities) == cMosaic.LCONE_ID);
    idxMconeCenter = find(theCenterMajorityConeType(mRGCindicesToVisualizeSTFsAcrossChromaticities) == cMosaic.MCONE_ID);

    % Visualize the locations of cells with a surround cone mix in the targetRangeForSurroundConeMix
    pdfFileNameLcenter = sprintf('SurroundConeMixOutlinedLcenterMRGClocationsWithinTargetRange_eccDegs_%2.1f_%2.1f.pdf', ...
        theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(1));
    pdfFileNameMcenter = sprintf('SurroundConeMixOutlinedMcenterMRGClocationsWithinTargetRange_eccDegs_%2.1f_%2.1f.pdf', ...
        theComputeReadyMRGCmosaic.eccentricityDegs(1), theComputeReadyMRGCmosaic.eccentricityDegs(1));


    hFigLcenter = figure(556); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theLcenterAxes = MSreadyPlot.generateAxes(hFigLcenter,ff, 'figPosition', [268 127]);

    hFigMcenter = figure(557); clf;
    theMcenterAxes = MSreadyPlot.generateAxes(hFigMcenter,ff, 'figPosition', [969 127]);

    plotTitleLcenter = sprintf('L-center mRGCs with surround cone mix in [%2.2f-%2.2f]',targetRangeForSurroundConeMix(1), targetRangeForSurroundConeMix(2));
    plotTitleMcenter = sprintf('M-center mRGCs with surround cone mix in [%2.2f-%2.2f]',targetRangeForSurroundConeMix(1), targetRangeForSurroundConeMix(2));
    MSreadyPlot.renderIdentifiedMRGClocations(hFigLcenter, hFigMcenter, theLcenterAxes{1,1}, theMcenterAxes{1,1}, ...
            theComputeReadyMRGCmosaic, ...
            mRGCindicesToVisualizeSTFsAcrossChromaticities(idxLconeCenter), ...
            mRGCindicesToVisualizeSTFsAcrossChromaticities(idxMconeCenter), ...
            plotTitleLcenter, plotTitleMcenter, ff);


    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileNameLcenter);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigLcenter, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileNameLcenter);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

    pdfFileNameUnscaled = fullfile(rawFiguresRoot, pdfFileNameMcenter);
    NicePlot.exportFigToPDF(pdfFileNameUnscaled, hFigMcenter, 300);

    if (exportScaledFigureVersionForManuscript)
        scaleFactor = 0.24;
        pdfFileNameScaled = fullfile(scaledFiguresRoot, pdfFileNameMcenter);
        MosaicPoolingOptimizer.exportScaledFigure(pdfFileNameUnscaled, pdfFileNameScaled, scaleFactor);
    end

end


function [surroundConeMixForAllCells, theCenterMajorityConeTypeForAllCells, centerConeMixForAllCells] = ...
    surroundConeMixForAllCellsInMosaic(theMRGCmosaic, performSurroundAnalysisForConesExclusiveToTheSurround)

    surroundConeMixForAllCells = zeros(1, theMRGCmosaic.rgcsNum);
    centerConeMixForAllCells = zeros(1, theMRGCmosaic.rgcsNum);
    theCenterMajorityConeTypeForAllCells = zeros(1, theMRGCmosaic.rgcsNum);

    parfor theRGCindex = 1:theMRGCmosaic.rgcsNum
        [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix, centerConeMix] = MosaicPoolingOptimizer.analyzeCenterSurroundConeMix(...
            theMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround);

        surroundConeMixForAllCells(theRGCindex) = surroundConeMix;
        centerConeMixForAllCells(theRGCindex) = centerConeMix;
        theCenterMajorityConeTypeForAllCells(theRGCindex) = theCenterMajorityConeType;
    end

end
