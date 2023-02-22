function visualizeSamplingGrids(obj)
   hFig = figure(2); clf;
   set(hFig, 'Color', [1 1 1]);

   conesNumPooledByTheRFcenters = unique(obj.conesNumPooledByTheRFcenterGrid);

   for iConesNumPooled = 1:numel(conesNumPooledByTheRFcenters)
        % Determine spatial grid coords for this # of center cones
        conesNumPooled = conesNumPooledByTheRFcenters(iConesNumPooled);
        gridNodesList = find(obj.conesNumPooledByTheRFcenterGrid == conesNumPooled);

        theLconeMajoritySpatialGrid = obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithLconeMajorityCenter(gridNodesList),:);
        theMconeMajoritySpatialGrid = obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithMconeMajorityCenter(gridNodesList),:);

        plotSamplingGridOnTopOfMosaic = false;
        ax = subplot(2,numel(conesNumPooledByTheRFcenters), iConesNumPooled);
        plotMultiFocalSpatialSamplingGrid(obj, ax, ...
            theLconeMajoritySpatialGrid, ...
            gridNodesList, ...
            plotSamplingGridOnTopOfMosaic, ...
            [1 0 0], ...
            sprintf('spatial sampling grid for %d L-center cones', conesNumPooled));

        ax = subplot(2,numel(conesNumPooledByTheRFcenters), iConesNumPooled+numel(conesNumPooledByTheRFcenters));
        plotMultiFocalSpatialSamplingGrid(obj, ax, ...
            theMconeMajoritySpatialGrid, ...
            gridNodesList, ...
            plotSamplingGridOnTopOfMosaic, ...
            [0 0.5 0], ...
            sprintf('spatial sampling grid for %d M-center cones', conesNumPooled));
    end

end

function plotMultiFocalSpatialSamplingGrid(obj, ax, spatialSamplingGrid, gridNodesList, plotSamplingGridOnTopOfMosaic, color, plotTitle)

    if (plotSamplingGridOnTopOfMosaic)
        obj.theRGCMosaic.visualize(...
            'axesHandle', ax, ...
            'samplingGrid', spatialSamplingGrid, ...
            'samplingGridOutlineColor', [1 1 0],...
            'samplingGridFillColor', [0 0 1],...
            'plotTitle', plotTitle);
    else
        plotMultiFocalSpatialSamplingGridIndices(obj, ax, spatialSamplingGrid, gridNodesList, color, plotTitle);
    end

    
end

function plotMultiFocalSpatialSamplingGridIndices(obj, ax, spatialSamplingGrid, gridNodesList, fontColor, plotTitle)
    centerDegs = obj.theRGCMosaic.eccentricityDegs;
    sizeDegs = obj.theRGCMosaic.sizeDegs;

    xx = centerDegs(1) + sizeDegs(1)/2*[-1 -1 1  1 -1];
    yy = centerDegs(2) + sizeDegs(2)/2*[-1  1 1 -1 -1];
    plot(ax,xx,yy,'k-', 'LineWidth', 1);
    hold(ax, 'on');
    for iNode = 1:numel(gridNodesList)
        text(ax, spatialSamplingGrid(iNode,1), spatialSamplingGrid(iNode,2), sprintf('%d', gridNodesList(iNode)), ...
            'FontSize', 14, 'Color', fontColor);
    end
    set(ax, 'FontSize', 16);
    axis(ax, 'equal')
    title(ax, plotTitle);
    box(ax, 'off');
    grid(ax, 'on');
    set(ax, 'XLim', centerDegs(1) + sizeDegs(1)*0.5*[-1 1] + [-0.1 0.1], ...
            'YLim', centerDegs(2) + sizeDegs(2)*0.5*[-1 1] + [-0.1 0.1]);
 
    drawnow;
end
