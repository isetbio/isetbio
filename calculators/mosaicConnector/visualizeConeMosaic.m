function visualizeConeMosaic(conePositionsMicrons, coneTypes, roi, plotlabOBJ)
    LconeIndices = find(coneTypes == 2);
    MconeIndices = find(coneTypes == 3);
    SconeIndices = find(coneTypes == 4);
    hFig = figure(); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.05);
    scatter(theAxesGrid{1,1}, conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 'r'); hold on
    scatter(theAxesGrid{1,1}, conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 'g');
    scatter(theAxesGrid{1,1}, conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 'b');
    LconePercent = 100*numel(LconeIndices)/(size(conePositionsMicrons,1));
    MconePercent = 100*numel(MconeIndices)/(size(conePositionsMicrons,1));
    SconePercent = 100*numel(SconeIndices)/(size(conePositionsMicrons,1));
    title(theAxesGrid{1,1}, sprintf('%2.2f%% (L), %2.2f%% (M), %2.2f%% (S)', LconePercent, MconePercent, SconePercent));
    
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);

    xLims = [xAxis(1) xAxis(end)] + roi.margin*[1,-1];
    yLims = [yAxis(1) yAxis(end)] + roi.margin*[1,-1];
    set(theAxesGrid{1,1}, 'CLim', [0 1], 'XLim', xLims, 'YLim', yLims);
    
    micronsPerDegree = 300;
    fName = sprintf('ConeMosaic_x=%2.2f_y=%2.2fdegs', roi.center(1)/micronsPerDegree, roi.center(2)/micronsPerDegree);
    plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(pwd(), 'exports'));
    
end
