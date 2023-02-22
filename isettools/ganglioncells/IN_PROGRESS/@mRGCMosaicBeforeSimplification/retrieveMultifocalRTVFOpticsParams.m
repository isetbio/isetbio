function retrieveMultifocalRTVFOpticsParams(obj, sourceMidgetRGCMosaic)
  
    % The number of center cones in the multifocal RTVF object
    centerConesNumExamined = sort(unique(sourceMidgetRGCMosaic.theConesNumPooledByTheRFcenterGrid), 'ascend');

    % The multifocalRTVFstruct
    obj.multifocalRTVFgrids = struct(...
        'conesNumPooledByTheRFcenterGrid', sourceMidgetRGCMosaic.theConesNumPooledByTheRFcenterGrid, ...
        'samplingPositionGrid', sourceMidgetRGCMosaic.theSamplingPositionGrid);
    obj.multifocalRTVFopticsParams = cell(1,numel(sourceMidgetRGCMosaic.theRetinaToVisualFieldTransformerOBJList));
    
    for iCone = 1:numel(centerConesNumExamined)
        iidx = find(sourceMidgetRGCMosaic.theConesNumPooledByTheRFcenterGrid == centerConesNumExamined(iCone));
        for idx = 1:numel(iidx)
            iObj = iidx(idx);
            if (obj.beVerbose)
                theLabel = sprintf('%d cones @ %+2.2f %+2.2f', ...
                    sourceMidgetRGCMosaic.theConesNumPooledByTheRFcenterGrid(iObj), ...
                    sourceMidgetRGCMosaic.theSamplingPositionGrid(iObj,1), ...
                    sourceMidgetRGCMosaic.theSamplingPositionGrid(iObj,2));
                fprintf('\nRTVF{%02d} was fitted for %s', iObj, theLabel);
            end

            theRTVFTobj = sourceMidgetRGCMosaic.theRetinaToVisualFieldTransformerOBJList{iObj};
            obj.multifocalRTVFopticsParams{iObj} = theRTVFTobj.opticsParams;
        end
    end

end

function visualizeMosaicBorders(obj, sourceMidgetRGCMosaic, theROI)
    
    sourceConesXYmin = min(obj.inputConeMosaic.coneRFpositionsDegs, [], 1);
    sourceConesXYmax = max(obj.inputConeMosaic.coneRFpositionsDegs, [], 1);
    
    theRGCXYmin = min(obj.rgcRFpositionsDegs, [],1);
    theRGCXYmax = max(obj.rgcRFpositionsDegs, [],1);

    sourceRGCXYmin = min(sourceMidgetRGCMosaic.rgcRFpositionsDegs, [],1);
    sourceRGCXYmax = max(sourceMidgetRGCMosaic.rgcRFpositionsDegs, [],1);


    hFig = figure(); clf;
    set(hFig, 'Position', [10 10 1500 800], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.08 0.5 0.90]);
    theROI.visualize(...
       'figureHandle', hFig, ...
       'axesHandle', ax, ...
       'fillColor', [0.8 0.8 0.2]);
    hold(ax, 'on');
    pCones = plotRect(ax, sourceConesXYmin(1), sourceConesXYmax(1), sourceConesXYmin(2), sourceConesXYmax(2), [0.5 0.5 0.5]);
    pSourceRGC = plotRect(ax, sourceRGCXYmin(1), sourceRGCXYmax(1), sourceRGCXYmin(2), sourceRGCXYmax(2), [0 0 1]);
    pRGC = plotRect(ax, theRGCXYmin(1), theRGCXYmax(1), theRGCXYmin(2), theRGCXYmax(2), [1 0 0]);
    legend(ax, [pCones pSourceRGC pRGC], {'input cone mosaic', 'source RGC mosaic', 'RGC mosaic'}, ...
        'Location', 'NorthOutside', 'Orientation', 'horizontal');
    axis(ax, 'equal');
    set(ax, 'XLim', [sourceConesXYmin(1)-0.2 sourceConesXYmax(1)+0.2], 'YLim', [sourceConesXYmin(2)-0.2 sourceConesXYmax(2)+0.2]);
    set(ax, 'XTick', -20:0.5:20, 'YTick', -20:0.5:20);
    
    grid(ax, 'on');
    xtickangle(ax, 0);
    set(ax, 'FontSize', 16, 'TickDir', 'both');
    xlabel(ax, 'x (degs)');
    ylabel(ax, 'y (degs)');

    
    ax = subplot('Position', [0.58 0.08 0.4 0.90]);
    pRGC = plotRect(ax, theRGCXYmin(1), theRGCXYmax(1), theRGCXYmin(2), theRGCXYmax(2), [1 0 0]);
    hold(ax, 'on');
    plot(obj.rgcRFpositionsDegs(:,1), obj.rgcRFpositionsDegs(:,2), 'ko');
    axis(ax, 'equal');
    set(ax, 'XLim', [theRGCXYmin(1)-0.2 theRGCXYmax(1)+0.2], 'YLim', [theRGCXYmin(2)-0.2 theRGCXYmax(2)+0.2]);
    set(ax, 'XTick', -20:0.2:20, 'YTick', -20:0.2:20);
    grid(ax, 'on');
    xtickangle(ax, 0);
    set(ax, 'FontSize', 16, 'TickDir', 'both');
    xlabel(ax, 'x (degs)');
    ylabel(ax, 'y (degs)');
    title('cropped mRGCmosaic')
end

function p = plotRect(ax, xMin, xMax, yMin, yMax, color)
    x = [xMin xMin xMax xMax xMin];
    y = [yMin yMax yMax yMin yMin];
    plot(ax, x, y, '-', 'Color', [0 0 0], 'LineWidth', 2);
    hold(ax, 'on')
    p = plot(ax, x, y, '--', 'Color', color, 'LineWidth', 2.0);
end


