function visualizeLattice(rfPositions)
    hFig = figure(1);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.04, ...
            'topMargin', 0.05);
    theAxesHandle = theAxesGrid{1,1};
    scatter(theAxesHandle, rfPositions(:,1), rfPositions(:,2), 2);
    axis(theAxesHandle, 'equal');
    xyRange = max(abs(rfPositions(:)))*[-1 1];
    set(theAxesHandle , 'XLim', xyRange, 'YLim', xyRange);
    plotlab.offsetAxes(theAxesHandle, 'offsetPercent', 0.015);
end
    
