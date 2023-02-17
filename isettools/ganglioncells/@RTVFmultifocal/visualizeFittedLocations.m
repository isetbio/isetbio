function visualizeFittedLocations(mosaicDirectory, figNo, theMidgetRGCmosaic, ...
    theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid)

    % Find how many different #of center cones are fitted
    conesNumPooled = unique(theConesNumPooledByTheRFcenterGrid);
    xCoords = squeeze(theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs(:,1));
    yCoords = squeeze(theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs(:,2));
    xLimsDegs = [min(xCoords) max(xCoords)];
    yLimsDegs = [min(yCoords) max(yCoords)];

    for iConesNumPooled = 1:numel(conesNumPooled)

        theConesNumPooled = conesNumPooled(iConesNumPooled);
        eccentricitySamplingGrid = RTVFmultifocal.subGridSpatialCoordsForConesNumPooled(theConesNumPooled, ...
            theConesNumPooledByTheRFcenterGrid, theOpticsPositionGrid);
        idx = find(sum(theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix,1) == theConesNumPooled);

        hFig = figure(figNo+iConesNumPooled); clf;
        ff = MSreadyPlot.figureFormat('1x1 large');
        set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
        
        % Left panel: RGC positions and RTVF sampling points (red)
        ax = subplot('Position', ff.subplotPosVectors(1, 1).v);
        MSreadyPlot.renderRFpositions(ax, theMidgetRGCmosaic.rgcRFpositionsDegs(idx,:), ...
            xLimsDegs, yLimsDegs, sprintf('%d-cone RF centers', theConesNumPooled), ff);

        hold(ax, 'on');
        plot(ax, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), '+', ...
           'Color', [1 0.5 0.5], 'LineWidth', ff.lineWidth*2, 'MarkerSize', ff.markerSize);
        plot(ax, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), '+', ...
           'Color', [1 0 0], 'LineWidth', ff.lineWidth, 'MarkerSize', ff.markerSize);

        drawnow

        pdfFilename = sprintf('RTVFmultifocal_%dCenterCone_samplingGrid.pdf', theConesNumPooled);
        pdfFilename = fullfile(mosaicDirectory,pdfFilename );
        NicePlot.exportFigToPDF(pdfFilename , hFig, 300);
        
    end
end
