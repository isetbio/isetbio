function peekIntoRTVFobj(theRTVFTobj, iRTVobjIndex, ...
    theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, figNo)

    iRTVobjIndexDoesNotCorrespondToFullList = false;
    if (iRTVobjIndex < 0)
        iRTVobjIndex = -iRTVobjIndex;
        iRTVobjIndexDoesNotCorrespondToFullList = true;
    end

    % Target and achieved ratios
    targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
    targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
    theRTVFTobj.rfComputeStruct.theSTF
    fittedRsRcRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedRsRcRatio;
    fittedSCintSensRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedSCIntSensRatio;
    targetKsKcRatio = targetSCintSensRatio / (targetRsRcRatio^2);
    fittedKsKcRatio = fittedSCintSensRatio / (fittedRsRcRatio^2);
    

    fprintf('RTVFobj at position (degs): %2.2f %2.2f with %d center cones\n', ...
        theSamplingPositionGrid(iRTVobjIndex,1), theSamplingPositionGrid(iRTVobjIndex,2), theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));
    fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
    fprintf('Target S/C int. sens. ratio: %2.3f, achieved: %2.3f\n', targetSCintSensRatio, fittedSCintSensRatio);


    figureName = sprintf('RTVF obj located at position (degs): %2.2f, %2.2f with %d center cones', ...
                         theSamplingPositionGrid(iRTVobjIndex,1), ...
                         theSamplingPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));

    RetinaToVisualFieldTransformer.performCronerKaplanSimulation(...
                theRTVFTobj.rfComputeStruct.theFittedVisualRFMap, ...
                theRTVFTobj.rfComputeStruct.theRetinalRFcenterConeMap, ...
                theRTVFTobj.rfComputeStruct.theRetinalRFsurroundConeMap, ...
                theRTVFTobj.rfComputeStruct.retinalConePoolingParams.finalValues, ...
                theRTVFTobj.rfComputeStruct.retinalConePoolingParams, ...
                theRTVFTobj.rfComputeStruct.theRotationOfTheFittedVisualRF, ...
                theRTVFTobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1), ...
                theRTVFTobj.rfComputeStruct.modelConstants.targetSTFmatchMode, ...
                theRTVFTobj.rfComputeStruct.theSTF.target, ...
                theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio, ...
                theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio, ...
                theRTVFTobj.rfComputeStruct.modelConstants.visualRFcenterCharacteristicRadiusDegs, ...
                theRTVFTobj.multiStartsNumDoGFit, ...
                figNo, figureName);

    ax = subplot(2,4,[7 8]);
    title(ax, sprintf('retinal pooling params for RFs with %d center cone(s) around (%2.2f,%2.2f) degs', ...
        theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), theSamplingPositionGrid(iRTVobjIndex,1), theSamplingPositionGrid(iRTVobjIndex,2)));

    % Plot correspondence between target and achieved ratios
    ax = subplot(2,4,4);
    yyaxis(ax, 'left');
    b1 = bar(ax, [2], 100*((fittedRsRcRatio/targetRsRcRatio)-1), 'BaseValue', 0.0);
    hold(ax, 'on');
    b2 = bar(ax, [4], 100*((fittedSCintSensRatio/targetSCintSensRatio)-1), 'BaseValue', 0.0);

    b1(1).FaceColor = 'flat';
    b1(1).CData(1,:) = [1 0 0];

    b2(1).FaceColor = 'flat';
    b2(1).CData(1,:) = [1 0 0];

    axis(ax, 'square');
    grid(ax, 'on');
    set(ax, 'YLim', [-1 1], 'YColor', [1 0 0], 'FontSize', 16, ...
        'YLim', [-100 100], 'YTick', -100:20:100);
    ylabel(ax,'(fitted-target)/target');

    yyaxis(ax, 'right')
    b = bar(ax, [6], 100*((fittedKsKcRatio/targetKsKcRatio)-1), 'BaseValue', 0.0);
    b(1).FaceColor = 'flat';
    b(1).CData(1,:) = [0 0 1];
    axis(ax, 'square');
    grid(ax, 'on');
    set(ax, 'XLim', [1 7],  'YColor', [0 0 1], 'XTick', [2 4 6], 'XTickLabel', {'Rs/Rc','S/C sens', 'Ks/Kc'}, ...
        'YLim', [-100 100], 'FontSize', 16, 'YTick', [-100:20:100]);
    ylabel(ax,'(fitted-target)/target');
    
    drawnow;

    hFig = figure(figNo);
    if (iRTVobjIndexDoesNotCorrespondToFullList)
        pdfFileName = sprintf('RTVFobj_X_at_%2.2f_%2.2fdegs_with%d_centerCones.pdf', ...
                         theSamplingPositionGrid(iRTVobjIndex,1), ...
                         theSamplingPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));
    else
        pdfFileName = sprintf('RTVFobj_%d_at_%2.2f_%2.2fdegs_with_%d_centerCones.pdf', ...
                         iRTVobjIndex, ...
                         theSamplingPositionGrid(iRTVobjIndex,1), ...
                         theSamplingPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));
    end
    
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

end