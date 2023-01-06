function peekIntoRTVFobj(theRTVFTobj, iRTVobjIndex, theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, figNo)
    % Target and achieved ratios
    targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
    targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
    fittedRsRcRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedRsRcRatio;
    fittedSCintSensRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedSCIntSensRatio;

    fprintf('RTVFobj at position (degs): %2.2f %2.2f with %d center cones\n', ...
        theSamplingPositionGrid(iRTVobjIndex,1), theSamplingPositionGrid(iRTVobjIndex,2), theConesNumPooledByTheRFcenterGrid(iRTVobjIndex));
    fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
    fprintf('Target S/C int. sens. ratio: %2.3f, achieved: %2.3f\n', targetSCintSensRatio, fittedSCintSensRatio);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 900 350], ...
        'Name', sprintf('RTVF obj #%d  located at position (degs): %2.2f, %2.2f with %d center cones', ...
                         iRTVobjIndex, ...
                         theSamplingPositionGrid(iRTVobjIndex,1), ...
                         theSamplingPositionGrid(iRTVobjIndex,2), ...
                         theConesNumPooledByTheRFcenterGrid(iRTVobjIndex)));
    ax = subplot(1,2,1);
    XLims = [1 15]; YLims = XLims;
    XTicks = 0:2:20; YTicks = XTicks;
    plotTargetAndAchievedParam(ax, targetRsRcRatio, fittedRsRcRatio, XLims, YLims, XTicks, YTicks, 'Rs/Rc ratio');

    ax = subplot(1,2,2);
    XLims = [0 1]; YLims = XLims;
    XTicks = 0:0.2:1.0; YTicks = XTicks;
    plotTargetAndAchievedParam(ax, targetSCintSensRatio, fittedSCintSensRatio, XLims, YLims, XTicks, YTicks, 'int. sens. S/C ratio');
    drawnow;
end

function plotTargetAndAchievedParam(ax, targetVal, fittedVal, XLims, YLims, XTicks, YTicks, titleString)
    
    plot(ax, XLims, YLims, 'k-', 'LineWidth', 1.0); hold(ax,'on');
    plot(ax, targetVal, fittedVal, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]);
    axis(ax, 'square');
    grid(ax, 'on');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', XTicks, 'YTick', YTicks, 'FontSize', 16);
    title(ax, titleString);
    xlabel(ax,'target');
    ylabel(ax,'achieved');
end