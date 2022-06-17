function visualizeRFoverlap2D(obj, iRGC, nearbyRGCindex, ...
    rfProfile2DmainRGC, rfProfile2DnearbyRGC, ...
    xSupport, ySupport, xTicks, yTicks, axRFOverlap2D)

    zLevels = [0.07 0.999];

    cmapMain = brewermap(1024, 'greys');
    alphaMain = 0.95;
    contourLineColorMain = [0 0 0];

    cmapNearby = brewermap(1024, 'reds');
    alphaNearby = 0.4;
    contourLineColorNearby = [1 0 0];
    hold(axRFOverlap2D, 'on');

    zData1 = (rfProfile2DmainRGC/max(rfProfile2DmainRGC(:))).^0.5;
    cMosaic.semiTransparentContourPlot(axRFOverlap2D, xSupport, ySupport, ...
        zData1, zLevels, cmapMain, alphaMain, contourLineColorMain, ...
        'lineWidth', 2.0);
    
    zData2 = (rfProfile2DnearbyRGC/max(rfProfile2DnearbyRGC(:))).^0.5;
    cMosaic.semiTransparentContourPlot(axRFOverlap2D, xSupport, ySupport, ...
              zData2, zLevels, cmapNearby, alphaNearby, contourLineColorNearby, ...
              'lineWidth', 2.0);

    % Compute the overlap coefficient
    overlapCoeff = RGCconnector.overlap(zData1(:), zData2(:));
    
    % Finalize plot
    axis(axRFOverlap2D,'equal');
    set(axRFOverlap2D, 'XLim', [xSupport(1) xSupport(end)], 'YLim', [ySupport(1) ySupport(end)], 'FontSize', 15);
    set(axRFOverlap2D, 'XTick', xTicks, 'YTick', yTicks);
    box(axRFOverlap2D, 'on'); grid(axRFOverlap2D, 'on');
    title(axRFOverlap2D, sprintf('Rc/RGCsep: %1.2f, 2D-overlap: %2.0f%%', ...
        obj.wiringParams.RcToRGCseparationRatio, 100*overlapCoeff));
    xlabel(axRFOverlap2D, 'space, x (microns)');
    ylabel(axRFOverlap2D, 'space, y (microns)');
end