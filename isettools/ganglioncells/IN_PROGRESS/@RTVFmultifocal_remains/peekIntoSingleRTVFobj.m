function peekIntoSingleRTVFobj(theRTVFTobj, iRTVobjIndex, ...
    theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, figNo)

    iRTVobjIndexDoesNotCorrespondToFullList = false;
    if (iRTVobjIndex < 0)
        iRTVobjIndex = -iRTVobjIndex;
        iRTVobjIndexDoesNotCorrespondToFullList = true;
    end

    peekIntoConeSpecificRTVFobj(...
        theRTVFTobj, ...
        cMosaic.LCONE_ID, ...
        theRTVFTobj.LconeRFcomputeStruct, ...
        iRTVobjIndexDoesNotCorrespondToFullList, ...
        iRTVobjIndex, ...
        theOpticsPositionGrid, ...
        theConesNumPooledByTheRFcenterGrid, ...
        figNo+1);


    peekIntoConeSpecificRTVFobj(...
        theRTVFTobj, ...
        cMosaic.MCONE_ID, ...
        theRTVFTobj.MconeRFcomputeStruct, ...
        iRTVobjIndexDoesNotCorrespondToFullList, ...
        iRTVobjIndex, ...
        theOpticsPositionGrid, ...
        theConesNumPooledByTheRFcenterGrid, ...
        figNo+2);


end

function peekIntoConeSpecificRTVFobj( ...
    theRTVFTobj, ...
    theConeType, ...
    theConeSpecificRFcomputeStruct, ...
    iRTVobjIndexDoesNotCorrespondToFullList, ...
    iRTVobjIndex, ...
    theOpticsPositionGrid, ...
    theConesNumPooledByTheRFcenterGrid, ...
    figNo)


    % Compute the visual RF map
    [theVisualRFmap, theRetinalRFcenterConeMap, ...
     theRetinalRFsurroundConeMap] = theRTVFTobj.visualRFfromRetinalConePooling(...
                theConeSpecificRFcomputeStruct.modelConstants, ...
                theConeSpecificRFcomputeStruct.retinalConePoolingParams.finalValues);
  
    if (strcmp(theRTVFTobj.stfComputeMethod, RTVF.modeledSTFcomputeMethod))
        % Compute the STF data from the visual RF map
        theSTFdata = theRTVFTobj.visualSTFfromCronerKaplanAnalysisOfVisualRF(theVisualRFmap, false);
        
    else
        theSTFdata = theConeSpecificRFcomputeStruct.theFinalSTFdata;
    end

    fittedRsRcRatio = theSTFdata.fittedDoGModelRsRcRatio
        fittedSCintSensRatio = theSTFdata.fittedDoGModelSCIntSensRatio

    % Target and achieved ratios
    targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
    targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

    
    targetKsKcRatio = targetSCintSensRatio / (targetRsRcRatio^2);
    fittedKsKcRatio = fittedSCintSensRatio / (fittedRsRcRatio^2);
    

    fprintf('RTVFobj at position (degs): %2.2f %2.2f with %d center cones\n', ...
        theOpticsPositionGrid(iRTVobjIndex,1), theOpticsPositionGrid(iRTVobjIndex,2), theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));
    fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
    fprintf('Target S/C int. sens. ratio: %2.3f, achieved: %2.3f\n', targetSCintSensRatio, fittedSCintSensRatio);


    switch (theConeType)
        case cMosaic.LCONE_ID
            conesString = 'L-cone(s)';
            
        case cMosaic.MCONE_ID
            conesString = 'M-cone(s)';
           
        otherwise
            error('Not an L- or an M-cone');
    end

    figureName = sprintf('RTVF obj %d located at position (degs): %2.2f, %2.2f with %d center %s', ...
                         iRTVobjIndex, ...
                         theOpticsPositionGrid(iRTVobjIndex,1), ...
                         theOpticsPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
                         conesString);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1000 800], 'Name', figureName, 'Color', [1 1 1]);

    ax = subplot(2,3,1);
    xSupportDegs = theConeSpecificRFcomputeStruct.modelConstants.spatialSupportDegs(:,1);
    ySupportDegs = theConeSpecificRFcomputeStruct.modelConstants.spatialSupportDegs(:,2);
    imagesc(xSupportDegs, ySupportDegs, theRetinalRFcenterConeMap)
    
    ax = subplot(2,3,2);
    imagesc(xSupportDegs, ySupportDegs, theRetinalRFsurroundConeMap)
 
    
    centerProfileX = sum(theRetinalRFcenterConeMap, 1);
    centerProfileY = sum(theRetinalRFcenterConeMap, 2);
    surroundProfileX = sum(theRetinalRFsurroundConeMap, 1);
    surroundProfileY = sum(theRetinalRFsurroundConeMap, 2);
    

    maxAll = max([max(centerProfileX) max(centerProfileY ) max(surroundProfileX ) max(surroundProfileY)]);
    maxSurround = max([max(surroundProfileX ) max(surroundProfileY)]);

    ax = subplot(2,3,3);
    plot(xSupportDegs, centerProfileX, 'r-', 'LineWidth', 1.5);
    hold on
    plot(xSupportDegs, -surroundProfileX, 'b-', 'LineWidth', 1.5);
    hold on
    plot(xSupportDegs, centerProfileX-surroundProfileX, 'k--', 'LineWidth', 1.0);
    set(gca, 'YLim', [-maxSurround maxAll]);
    

    ax = subplot(2,3,4);
    plot(theSTFdata.spatialFrequencySupport, theSTFdata.visualSTF, 'ks'); hold on
    plot(theSTFdata.spatialFrequencySupport, theSTFdata.fittedDoGModelToVisualSTF.compositeSTF, 'k-');
    plot(theSTFdata.spatialFrequencySupport, theSTFdata.fittedDoGModelToVisualSTF.centerSTF, 'r-');
    plot(theSTFdata.spatialFrequencySupport, theSTFdata.fittedDoGModelToVisualSTF.surroundSTF, 'b-');
    set(gca, 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100], 'XLim', [0.1 100]);
    title(ax, sprintf('Rs/Rc: %2.2f (target: %2.2f, mismatch: %2.1f%%), S/CintSens: %2.2f (target: %2.2f, mismatch:%2.1f%%)', ...
        fittedRsRcRatio, targetRsRcRatio, 100*((fittedRsRcRatio/targetRsRcRatio)-1), ...
        fittedSCintSensRatio, targetSCintSensRatio, 100*((fittedSCintSensRatio/targetSCintSensRatio)-1)));

    % Plot correspondence between target and achieved ratios
    ax = subplot(2,3,5);
    bar(ax, 1, 100*((fittedRsRcRatio/targetRsRcRatio)-1), 'BaseValue', 0.0);
    hold(ax, 'on');
    bar(ax, 2, 100*((fittedSCintSensRatio/targetSCintSensRatio)-1), 'BaseValue', 0.0);


    bar(ax, 3, 100*((fittedKsKcRatio/targetKsKcRatio)-1), 'BaseValue', 0.0);
    axis(ax, 'square');
    grid(ax, 'on');
    set(ax, 'XLim', [0 4],  'XTick', [1 2 3], 'XTickLabel', {'Rs/Rc','S/C sens', 'Ks/Kc'}, ...
        'YLim', [-100 100], 'FontSize', 16, 'YTick', -100:20:100);
    ylabel(ax,'(fitted-target)/target');
    
    ax = subplot(2,3,6);
    plot(ySupportDegs, centerProfileY, 'r-', 'LineWidth', 1.5);
    hold on
    plot(ySupportDegs, -surroundProfileY, 'b-', 'LineWidth', 1.5);
    plot(ySupportDegs, centerProfileY-surroundProfileY, 'k--', 'LineWidth', 1.0);
    hold on
    set(gca, 'YLim', [-maxSurround maxAll]);

    colormap(brewermap(1024, 'greys'));
    drawnow;
    
    if (iRTVobjIndexDoesNotCorrespondToFullList)
        pdfFileName = sprintf('RTVFobj_X_at_%2.2f_%2.2fdegs_with%d_center_%s.pdf', ...
                         theOpticsPositionGrid(iRTVobjIndex,1), ...
                         theOpticsPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
                         conesString);
    else
        pdfFileName = sprintf('RTVFobj_%d_at_%2.2f_%2.2fdegs_with_%d_center_%s.pdf', ...
                         iRTVobjIndex, ...
                         theOpticsPositionGrid(iRTVobjIndex,1), ...
                         theOpticsPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
                         conesString);
    end
    
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

end