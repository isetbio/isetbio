function hFig = plotMosaicProgression(obj)

    hFig = figure(1); 
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1050 1350]);
    
    rowsNum = 4;
    colsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.005, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.01);
       
    row = 1; col = 1; iteration = 0;
    plotMosaic(obj, subplotPosVectors, row, col, iteration);
    
    row = 1; col = 2; iteration = 1;
    plotMosaic(obj, subplotPosVectors, row, col, iteration);
    
    row = 3; col = 1; iteration = 100;
    plotMosaic(obj, subplotPosVectors, row, col, iteration);
    
    row = 3; col = 2; iteration = size(obj.latticeAdjustmentSteps,1);
    plotMosaic(obj, subplotPosVectors, row, col, iteration);
    
    NicePlot.exportFigToPDF('mosaicConstruction.pdf', hFig, 300);
    
    % Restore the final state of coneLocsHexGrid
    obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(end,:,:));
end

function plotMosaic(obj, subplotPosVectors, row, col, iteration)

    % Get cone locs at desired iteration
    if (iteration == 0)
        obj.coneLocsHexGrid = obj.initialLattice;
    else
        obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(iteration,:,:));
    end
    
    % Compute actual mosaic density
    [densityMapMosaic, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
    
    % Compute model mosaic density
    [densityMapModel, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
    
    
    % Positions in meters
    xRange = max(abs(obj.coneLocsHexGrid(:)))* [-1 1];
    displayedXrange = xRange(end)*0.76;
    displayedYrange = xRange(2)*0.23;
    
    % Center of focal, magnified view of the mosaic
    focusXpoint =  xRange(2)*0.17;
    focusYpoint = -xRange(2)*0.15;
    % Extent of the magnified view
    focusXrange = 2e-6*([-18 18]);
    focusYrange = 2e-6*([-11 11]);
    displayedXrange2 = focusXpoint + focusXrange;
    displayedYrange2 = focusYpoint + focusYrange;
    focusAreaOutline.x = [displayedXrange2(1) displayedXrange2(1) displayedXrange2(2) displayedXrange2(2) displayedXrange2(1)];
    focusAreaOutline.y = [displayedYrange2(1) displayedYrange2(2) displayedYrange2(2) displayedYrange2(1) displayedYrange2(1)];
    
    subplot('Position', subplotPosVectors(row,col).v);
    idx = find(...
            (obj.coneLocsHexGrid(:,1) > 1e-6) & ...
            (obj.coneLocsHexGrid(:,1) <= displayedXrange) & ...
            (abs(obj.coneLocsHexGrid(:,2)) <= displayedYrange)...
        );
    plot(obj.coneLocsHexGrid(idx,1), obj.coneLocsHexGrid(idx,2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    hold on; 
    
     % plot density contour
    if (iteration> 0)
        contourLevels = 1000*(25: 25: 240);
        plotContoursOverHalfField = false;
        if (plotContoursOverHalfField)
            idx = find(~(...
                (densityMapSupportX >= -1.5*1e-6) & ...
                (densityMapSupportX <= 2*displayedXrange) & ...
                (densityMapSupportY >= -1.5*1e-6-displayedYrange) & ...
                (densityMapSupportY >= displayedYrange)));
            densityMap(idx) = NaN;
        end
    
        contour(gca, densityMapSupportX, densityMapSupportY, densityMapMosaic, contourLevels, 'LineColor', [1.0 0.5 0.5], 'LineWidth', 3.0, 'ShowText', 'off');
        [cH, hH] = contour(gca, densityMapSupportX, densityMapSupportY, densityMapMosaic, contourLevels, 'LineColor', 'r', 'LineWidth', 1.5, 'ShowText', 'on', 'LabelSpacing',600);
        clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [1.0 0.2 0.2], 'Background', [1 1 1]);
   
        contour(gca, densityMapSupportX, densityMapSupportY, densityMapModel, contourLevels, 'LineColor', [0.5 0.5 1.0], 'LineWidth', 3.0, 'ShowText', 'off');
        [cH, hH] = contour(gca, densityMapSupportX, densityMapSupportY, densityMapModel, contourLevels, 'LineColor', 'b', 'LineWidth', 1.5, 'ShowText', 'on', 'LabelSpacing',300);
        clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [0.2 .2 1.0], 'Background', [1 1 1]);
    end
    
    %set(gca, 'CLim', [contourLevels(1) contourLevels(end)]);
    
    % Plot the focus area
    plot(focusAreaOutline.x, focusAreaOutline.y, 'k-', 'LineWidth', 1.5);
    hold off;
    
    axis 'equal';
    axis 'xy';
    set(gca,  'XLim', [-1.5*1e-6 displayedXrange], 'YLim', [-2*1e-6-displayedYrange displayedYrange], 'FontSize', 14, 'TickDir', 'both',  'LineWidth', 1.0);
    xTicks = 1e-6*(0:50:1000);
    if (col == 1)
        yTicks = 1e-6*(-1000:50:1000);
        yAxisColor = 'k';
    else
        yTicks = [];
        yAxisColor = 'none';
    end
    set(gca, 'YColor', yAxisColor, 'XTick', xTicks, 'XTickLabels', sprintf('%2.0f\n', xTicks/1e-6), 'YTick', yTicks, 'YTickLabels', sprintf('%2.0f\n', yTicks/1e-6));
    grid on; box off;
    title(sprintf('\niteration %d', iteration), 'FontSize', 18, 'FontWeight', 'bold');
    
    
    subplot('Position', subplotPosVectors(row+1,col).v);
%     idx = find(...
%         (obj.coneLocsHexGrid(:,1) > displayedXrange2(1)-6*1e-6) & ...
%         (obj.coneLocsHexGrid(:,1) < displayedXrange2(2)+6*1e-6) & ...
%         (obj.coneLocsHexGrid(:,2) > displayedYrange2(1)-6*1e-6) & ...
%         (obj.coneLocsHexGrid(:,2) < displayedYrange2(2)+6*1e-6)...
%         );
    
%    renderHexMesh(gca, obj.coneLocsHexGrid(idx,1), obj.coneLocsHexGrid(idx,2), [0.5 0.5 0.5], 'none', 0.8, 0.2, '-');
    hold on;
    iTheta = (0:10*360)/180*pi;
    coneApertureRadius = 0.98;
    xAperture = coneApertureRadius*cos(iTheta)*1e-6;
    yAperture = coneApertureRadius*sin(iTheta)*1e-6;
    
    idx = find(...
        (obj.coneLocsHexGrid(:,1) > displayedXrange2(1)-3*1e-6) & ...
        (obj.coneLocsHexGrid(:,1) < displayedXrange2(2)+3*1e-6) & ...
        (obj.coneLocsHexGrid(:,2) > displayedYrange2(1)-3*1e-6) & ...
        (obj.coneLocsHexGrid(:,2) < displayedYrange2(2)+3*1e-6)...
        );
    edgeColor = [0 0 0];
    faceColor = 0.85* [1 1 1];
    edgeAlpha = 1;
    faceAlpha = 1;
    for k = 1:numel(idx)
        x = squeeze(obj.coneLocsHexGrid(idx(k),1))+xAperture;
        y = squeeze(obj.coneLocsHexGrid(idx(k),2))+yAperture;
        patch(x, y, [0 0 0], 'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha, 'FaceAlpha', faceAlpha, 'FaceColor', faceColor, 'LineWidth', 1.5, 'LineStyle', '-', 'Parent', gca);
        plot(x, y, 'k-', 'LineWidth', 1.5);
    end
    
    if (iteration> 0)
        contour(gca, densityMapSupportX, densityMapSupportY, densityMapMosaic, contourLevels, 'LineColor', [1.0 0.5 0.5], 'LineWidth', 3.0, 'ShowText', 'off');
        [cH, hH] = contour(gca, densityMapSupportX, densityMapSupportY, densityMapMosaic, contourLevels, 'LineColor', 'r', 'LineWidth', 1.5, 'ShowText', 'on', 'LabelSpacing',300);
        clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [1.0 0.2 0.2], 'Background', [1 1 1]);
        contour(gca, densityMapSupportX, densityMapSupportY, densityMapModel, contourLevels, 'LineColor', [0.5 0.5 1.0], 'LineWidth', 3.0, 'ShowText', 'off');
        [cH, hH] = contour(gca, densityMapSupportX, densityMapSupportY, densityMapModel, contourLevels, 'LineColor', 'b', 'LineWidth', 1.5, 'ShowText', 'on', 'LabelSpacing',300);
        clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [0.2 .2 1.0], 'Background', [1 1 1]);
    end
    
    hold off
    axis 'equal';
    axis 'xy'
    box 'on'
    set(gca, 'XLim', [displayedXrange2(1) displayedXrange2(2)+1*1e-6], 'YLim', [displayedYrange2(1)-0.5*1e-6 displayedYrange2(2)+1.5*1e-6]);
   
    xTicks = 1e-6*(-500:10:500);
    if (col == 1)
        yTicks = 1e-6*(-200:10:200);    
    else
        yTicks = [];
        xTicks = [];
    end
    yAxisColor = 'k';
    set(gca, 'FontSize', 14, 'YColor', yAxisColor, 'XTick', xTicks, 'XTickLabels', sprintf('%2.0fum\n', xTicks/1e-6), 'YTick', yTicks, 'YTickLabels', sprintf('%2.0f\n', yTicks/1e-6));
    xlabel('microns');
    drawnow
    
    
end

function renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
    x = []; y = [];
    triangleConeIndices = delaunayn([xHex(:), yHex(:)]);
    for triangleIndex = 1:size(triangleConeIndices,1)
        coneIndices = triangleConeIndices(triangleIndex, :);
        xCoords = xHex(coneIndices);
        yCoords = yHex(coneIndices);
        for k = 1:numel(coneIndices)
            x = cat(2, x, xCoords);
            y = cat(2, y, yCoords);
        end
    end
    patch(x, y, [0 0 0], 'EdgeColor', meshEdgeColor, 'EdgeAlpha', meshEdgeAlpha, 'FaceAlpha', meshFaceAlpha, 'FaceColor', meshFaceColor, 'LineWidth', 1.5, 'LineStyle', lineStyle, 'Parent', axesHandle);
end
