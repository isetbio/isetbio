function visualizeVisualRFmapForTargetRGC(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType)

    % Find the target RGC to be visualized
    [targetCenterConesNumNotMatched, theVisualizedRGCindex] = theComputeReadyMRGCmosaic.indexOfRGCNearPosition( ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);
    if (targetCenterConesNumNotMatched)
        fprintf(2, 'Could not find an RGC with %d center cones\n', targetCenterConesNum);
        return;
    end

    % Load the optimall mapped visual RF maps for all cells
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps');

    % Figure format
    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x4 RF poster');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    tickSeparationArcMin = 6;
    spatialSupportRangeArcMin = tickSeparationArcMin*4;

    % Plot the visual RF map
    retinalRGCRFposDegs = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(theVisualizedRGCindex,:);
    mRGCMosaic.visualizeVisualRFmap(...
        optimallyMappedVisualRFmaps{theVisualizedRGCindex}, ...
        retinalRGCRFposDegs, ...
        theAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
        'withFigureFormat', ff);

end

