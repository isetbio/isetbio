function visualizeRGCmosaic(figNo, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, roi, stateLabel, plotlabOBJ)

    xOutline = cosd(0:10:360);
    yOutline = sind(0:10:360);
      
    hFig = figure(figNo); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.05);
    theAxesGrid = theAxesGrid{1,1};
    hold(theAxesGrid, 'on');
    
    % Sampling for contours
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    for iRGC = 1:rgcsNum
        
        xPts = RGCRFPositionsMicrons(iRGC,1)+0.5*RGCRFSpacingsMicrons(iRGC)*xOutline;
        yPts = RGCRFPositionsMicrons(iRGC,2)+0.5*RGCRFSpacingsMicrons(iRGC)*yOutline;
        patch(theAxesGrid,xPts, yPts,[0.8 0.8 0.8]); 
    end
    
    title(theAxesGrid,sprintf('RGC mosaic (%s)', stateLabel));
    
    xLims = [xAxis(1) xAxis(end)] + roi.margin*[1,-1];
    yLims = [yAxis(1) yAxis(end)] + roi.margin*[1,-1];
   
    set(theAxesGrid, 'CLim', [0 1], 'XLim', xLims, 'YLim', yLims);
    
    % Export the figure to the gallery directory in PNG format
    micronsPerDegree = 300;
    fName = sprintf('RGCMosaic(%s)_x=%2.2f_y=%2.2fdegs', stateLabel, roi.center(1)/micronsPerDegree, roi.center(2)/micronsPerDegree);
    plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(pwd(), 'exports'));
    
end

