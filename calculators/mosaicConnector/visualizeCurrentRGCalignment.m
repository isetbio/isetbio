function visualizeCurrentRGCalignment(rgcPMicrons, rgcPMicronsAligned, RGCRFSpacingMicrons, ...
                    xOutline, yOutline, coneIndicesWithinReach, conePositionsMicrons, coneSpacingsMicrons, desiredConesToRGCratios)
    hFig = figure(1); clf;
    theAxesGrid = plotlab.axesGrid(hFig, 'leftMargin', 0.03);

    xPts = rgcPMicrons(1)+0.5*RGCRFSpacingMicrons*xOutline;
    yPts = rgcPMicrons(2)+0.5*RGCRFSpacingMicrons*yOutline;
    plot(theAxesGrid{1,1}, rgcPMicrons(1)+0.5*RGCRFSpacingMicrons*xOutline, ...
         rgcPMicrons(2)+0.5*RGCRFSpacingMicrons*yOutline, 'r-');
    hold(theAxesGrid{1,1}, 'on')

    plot(theAxesGrid{1,1}, rgcPMicronsAligned(1)+0.5*RGCRFSpacingMicrons*xOutline, ...
         rgcPMicronsAligned(2)+0.5*RGCRFSpacingMicrons*yOutline, 'g-');

    for k = 1:numel(coneIndicesWithinReach)
        coneIndex = coneIndicesWithinReach(k);
        xPts = cat(2, xPts, conePositionsMicrons(coneIndex,1) + 0.5*coneSpacingsMicrons(coneIndex)*xOutline);
        yPts = cat(2, yPts, conePositionsMicrons(coneIndex,2) + 0.5*coneSpacingsMicrons(coneIndex)*yOutline);
        plot(theAxesGrid{1,1}, conePositionsMicrons(coneIndex,1) + 0.6*coneSpacingsMicrons(coneIndex)*xOutline, ...
             conePositionsMicrons(coneIndex,2) + 0.6*coneSpacingsMicrons(coneIndex)*yOutline,'b-');
    end

    hold(theAxesGrid{1,1}, 'off')
    xMin = min(xPts); xMax = max(xPts);
    yMin = min(yPts); yMax = max(yPts);

    xLim = [xMin xMax];
    yLim = [yMin yMax];

    set(theAxesGrid{1,1}, 'XLim', xLim, 'YLim', yLim);
    axis(theAxesGrid{1,1}, 'equal')
    title(theAxesGrid{1,1},sprintf('cone-to-RGC ratio: %2.2f', desiredConesToRGCratios(rgcIndex)));
end
