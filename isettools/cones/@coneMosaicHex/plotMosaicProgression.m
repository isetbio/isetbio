function hFig = plotMosaicProgression(obj, varargin)

    p = inputParser;
    p.addParameter('contourLevels', 1e3*[150 175 200 225 250], @isnumeric);
    p.addParameter('intermediateIterationsToDisplay', [10 100], @(x)(isnumeric(x)&&numel(x)==2));
    p.addParameter('displayedXrangeDegs', [], @isnumeric);
    p.addParameter('displayedYrangeDegs', [], @isnumeric);
    p.parse(varargin{:});
    
    % Get params
    contourLevels = p.Results.contourLevels;
    intermediateIterationsToDisplay = p.Results.intermediateIterationsToDisplay;
    
    micronsPerDegree = 300;
    displayedXrangeMicrons = p.Results.displayedXrangeDegs * micronsPerDegree;
    displayedYrangeMicrons = p.Results.displayedYrangeDegs * micronsPerDegree;
    
    hFig = figure(1);  clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1360 950]);
    
    rowsNum = 2;
    colsNum = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin',   0.06, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.02);
       
    row = 1; col = 1; iteration = 0; 
    labelContours = false; labelCones = false; backgroundColor = [1 1 1];
    subplotTitle = 'regular hexagonal lattice';
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, false, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle);
    
    row = 1; col = 2; iteration = 1; 
    labelContours = false; labelCones = false;
    subplotTitle = 'subsampling'; 
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, false, contourLevels,labelContours,  labelCones, backgroundColor, iteration, subplotTitle );
    
    row = 1; col = 3; iteration = intermediateIterationsToDisplay(1);  
    labelContours = false; labelCones = false;
    subplotTitle = sprintf('iteration #%d', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, true, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle );
    
    row = 2; col = 1; iteration = intermediateIterationsToDisplay(2);  
    labelContours = false; labelCones = false; 
    subplotTitle = sprintf('iteration #%d', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, true, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle );
    
    row = 2; col = 2; iteration = size(obj.latticeAdjustmentSteps,1); 
    labelContours = true; labelCones = false; 
    subplotTitle = sprintf('iteration #%d (converged)', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, true, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle );
    
    row = 2; col = 3; iteration = size(obj.latticeAdjustmentSteps,1); 
    labelContours = false; labelCones = true; contourLevels = [];
    subplotTitle = 'cone type assignment'; 
    plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, false, true, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle );
    
    
    % Restore the final state of coneLocsHexGrid
    obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(end,:,:));
end


function plotMosaic(obj, subplotPosVectors, row, col, displayedXrangeMicrons, displayedYrangeMicrons, plotHexMesh, plotContours, contourLevels, labelContours, labelCones, backgroundColor, iteration, subplotTitle)
    % Get cone locs at desired iteration
    if (iteration == 0)
        obj.coneLocsHexGrid = obj.initialLattice;
    else
        obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(iteration,:,:));
    end
    
    if (isempty(displayedXrangeMicrons))
        displayedXrangeMeters = obj.center(1)+ obj.width/2*[-1 1];
    else
        displayedXrangeMeters = displayedXrangeMicrons/2 * 1e-6 * [-1 1];
    end
    if (isempty(displayedXrangeMicrons))
        displayedYrangeMeters = obj.center(2)+ obj.height/2*[-1 1];
    else
        displayedYrangeMeters = displayedYrangeMicrons/2 * 1e-6 * [-1 1];
    end
    
    % Determine cones to be plotted
    coneApertureMicrons = obj.lambdaMin/1e6;
    idx = find(...
            (obj.coneLocsHexGrid(:,1)-0.5*coneApertureMicrons > displayedXrangeMeters(1)) & ...
            (obj.coneLocsHexGrid(:,1)+0*coneApertureMicrons < displayedXrangeMeters(2)) & ...
            (obj.coneLocsHexGrid(:,2)-0.5*coneApertureMicrons > displayedYrangeMeters(1)) & ...
            (obj.coneLocsHexGrid(:,2)+0*coneApertureMicrons < displayedYrangeMeters(2)) ...
        );
    
    ax = axes('Position', subplotPosVectors(row,col).v, 'units', 'pixels', 'Color', backgroundColor);
        
    if (labelCones)
        obj.visualizeGrid('axesHandle', ax, 'visualizedConeAperture', 'both');
    else   
        % Plot the cones
        edgeColor = 'none';
        faceColor = 0.7*[1 1 1];
        lineStyle = '-';
        iTheta = (0:60:300)/180*pi;
        coneApertureRadius = obj.lambdaMin/2;
        coneAperture.x = coneApertureRadius*cos(iTheta)*1e-6;
        coneAperture.y = coneApertureRadius*sin(iTheta)*1e-6;
        coneMosaicHex.renderPatchArray(ax, coneAperture, squeeze(obj.coneLocsHexGrid(idx,1)), squeeze(obj.coneLocsHexGrid(idx,2)), edgeColor, faceColor, lineStyle);
    end
    
    hold on
    
    % superimpose density contours
    if (plotContours > 0) && (~isempty(contourLevels))
        % Compute actual mosaic density
        [densityMapMosaic, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from mosaic');
    
        % Compute model mosaic density
        [densityMapModel, densityMapSupportX, densityMapSupportY] = obj.computeDensityMap('from model');
   
        [cH, hH] = contour(ax, densityMapSupportX, densityMapSupportY, densityMapModel, contourLevels, ...
                'LineColor', [0.1 .1 1.0], 'LineWidth', 3.0, 'ShowText', 'off');
        if (labelContours)
            clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [0.1 0.1 1.0], 'Background', [1 1 1], 'LabelSpacing', 300);
        end
            
        [cH, hH] = contour(ax, densityMapSupportX, densityMapSupportY, densityMapMosaic, contourLevels, ...
                'LineColor', [1 0.1 0.1], 'LineWidth', 3.0, 'ShowText', 'off');
        if (labelContours)
            clabel(cH,hH,'FontWeight','bold', 'FontName', 'Menlo', 'FontSize', 12, 'Color', [1 0.1 0.1], 'Background', [1 1 1], 'LabelSpacing', 100);
        end
    end
    
    if (plotHexMesh)
        coneMosaicHex.renderHexMesh(ax, obj.coneLocsHexGrid(idx,1), obj.coneLocsHexGrid(idx,2), [0.5 0.5 0.5], 'none', 0.8, 0.2, '-');
    end
    
    % Finalize plot
    hold off;
    axis 'equal';
    axis 'xy';
    xTicks = 1e-6*(-1000:50:1000);
    if (row == 1)
        xTickLabels = {};
    else
        xTickLabels = sprintf('%2.0f\n', xTicks/1e-6);
    end
    if (col == 1)
        yTicks = 1e-6*(-1000:50:1000);
        yAxisColor = 'k';
        ylabel('microns');
    else
        yTicks = [];
        yAxisColor = 'none';
    end
    yTickLabels = sprintf('%2.0f\n', yTicks/1e-6);
    set(gca,  ...
        'XLim', [displayedXrangeMeters(1)-2*1e-6 displayedXrangeMeters(2)+3*1e-6], ...
        'YLim', [displayedYrangeMeters(1)-2*1e-6 displayedYrangeMeters(2)+3*1e-6], ...
        'XTick', xTicks, 'XTickLabels', xTickLabels, ...
        'YTick', yTicks, 'YTickLabels', yTickLabels, ...
        'TickDir', 'both',  ...
        'FontSize', 18, ...
        'YColor', yAxisColor, ...
        'LineWidth', 1.0);
    grid off; box off;
    
    title(subplotTitle, 'FontSize', 18, 'FontWeight', 'bold');
    drawnow;
end