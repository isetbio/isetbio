function visualizeLatticeAndQuality(rfPositions, spacingDeviations, reTriangulationIsNeeded, iteration)
    hFig = figure(2); clf;
   
    theAxesGrid = plotlab.axesGrid(hFig, ...
        'rowsNum', 1, ...
        'colsNum', 2, ...
        'leftMargin', 0.04, ...
        'widthMargin', 0.07, ...
        'heightMargin', 0.07, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.06, ...
        'topMargin', 0.05);
    
    % The lattice on the left
    theAxesHandle = theAxesGrid{1,1};
    plot(theAxesHandle, rfPositions(:,1), rfPositions(:,2), 'k.'); 
    hold(theAxesHandle, 'on');
    [min(spacingDeviations) max(spacingDeviations)]
    idx = find(spacingDeviations < 0.5);
    plot(theAxesHandle, rfPositions(idx,1), rfPositions(idx,2), 'r.'); 
    idx = find(spacingDeviations < 0.3);
    plot(theAxesHandle, rfPositions(idx,1), rfPositions(idx,2), 'm.'); 
    hold(theAxesHandle, 'off');
    axis(theAxesHandle, 'equal');
    xyRange = max(abs(rfPositions(:)))*[-1 1];
    axis(theAxesHandle, 'square');
    set(theAxesHandle , 'XLim', xyRange, 'YLim', xyRange);
    plotlab.offsetAxes(theAxesHandle, 'offsetPercent', 0.015);   
    if (reTriangulationIsNeeded)
        title(theAxesHandle,'Re-triangularization in next iteration', 'Color', 'r');
    end
    
    % The quality on the right
    theAxesHandle = theAxesGrid{1,2};
    triangleIndices = delaunayn(rfPositions);
    [minQualityValue, histogramData] = computeHexLatticeQuality(rfPositions, triangleIndices);
    bar(theAxesHandle, histogramData.x, histogramData.y, 'FaceColor', [0.5 0.5 0.5]); hold on;
    plot(theAxesHandle, minQualityValue*[1 1], [0 max(histogramData.y)], 'r');
    grid(theAxesHandle, 'on')
    axis(theAxesHandle, 'square');
    xlabel(theAxesHandle, 'hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex');
    title(theAxesHandle,sprintf('iteration %d', iteration));
end
    
