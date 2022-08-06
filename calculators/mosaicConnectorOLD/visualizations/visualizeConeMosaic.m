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
    scatter(theAxesGrid{1,1}, conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 25, 'r'); hold on
    scatter(theAxesGrid{1,1}, conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 25, 'g');
    scatter(theAxesGrid{1,1}, conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 25, 'b');
    LconePercent = 100*numel(LconeIndices)/(size(conePositionsMicrons,1));
    MconePercent = 100*numel(MconeIndices)/(size(conePositionsMicrons,1));
    SconePercent = 100*numel(SconeIndices)/(size(conePositionsMicrons,1));
    title(theAxesGrid{1,1}, sprintf('%2.2f%% (L), %2.2f%% (M), %2.2f%% (S)', LconePercent, MconePercent, SconePercent));
    
    
    deltaX = min([roi.width/2 roi.height/2]);
    xAxis = 0 : deltaX: (roi.xo+roi.width/2);
    yAxis = 0 : deltaX: (roi.yo+roi.height/2);
    

    xAxisMicrons = 1000*WatsonRGCModel.rhoDegsToMMs(xAxis);
    yAxisMicrons = 1000*WatsonRGCModel.rhoDegsToMMs(yAxis);
    xAxisMicrons = [-fliplr(xAxisMicrons(2:end)) xAxisMicrons];
    yAxisMicrons = [-fliplr(yAxisMicrons(2:end)) yAxisMicrons];
    
    xLims = [xAxisMicrons(1) xAxisMicrons(end)];
    yLims = [yAxisMicrons(1) yAxisMicrons(end)];
    

    axis(theAxesGrid{1,1}, 'equal');
    set(theAxesGrid{1,1}, 'CLim', [0 1], 'XLim', xLims, 'YLim', yLims);

    fName = sprintf('ConeMosaic_x=%2.2f_y=%2.2fdegs', roi.xo, roi.yo);
   % plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(pwd(), 'exports'));
    
end
