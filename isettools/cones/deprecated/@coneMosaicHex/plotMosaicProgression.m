function hFig = plotMosaicProgression(obj, varargin)
% Plot the mosaic progression
%
% Syntax:
%   hFig = plotMosaicProgression(obj, [varargin])
%
% Description:
%    Plot the mosaic progression
%
% Inputs:
%    obj - the cone mosaic hex object
%
% Outputs:
%    hFig - The figure handle
%
% Optional key/value pairs:
%  
%  'contourLevels'                    - Vector. cone density levels at which to draw contours
%  'intermediateIterationsToDisplay'  - 2 element vector. 2 iterations at 
%                                       which to display the in-progress cone mosaic
%   'displayedXrangeDegs'             - visualzed x-axis range in degrees
%   'displayedYrangeDegs'             - visualzed y-axis range in degrees
%

    p = inputParser;
    p.addParameter('contourLevels', 1e3 * [150 175 200 225 250], ...
        @isnumeric);
    p.addParameter('intermediateIterationsToDisplay', [10 100], ...
        @(x)(isnumeric(x) && numel(x) == 2));
    p.addParameter('displayedXrangeDegs', [], @isnumeric);
    p.addParameter('displayedYrangeDegs', [], @isnumeric);
    p.parse(varargin{:});

    % Get params
    contourLevels = p.Results.contourLevels;
    intermediateIterationsToDisplay = ...
        p.Results.intermediateIterationsToDisplay;

    micronsPerDegree = 300;
    displayedXrangeMicrons = p.Results.displayedXrangeDegs * ...
        micronsPerDegree;
    displayedYrangeMicrons = p.Results.displayedYrangeDegs * ...
        micronsPerDegree;
    hFig = figure(1);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 620 1125]);

    rowsNum = 6;
    colsNum = 1;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.015, ...
           'widthMargin', 0.00, ...
           'leftMargin', 0.04, ...
           'rightMargin', -0.02, ...
           'bottomMargin', 0.05, ...
           'topMargin', 0.00);

    backgroundColor = 0.7*[1 1 1];

    col = 1;
    row = 1;
    iteration = 0;
    labelContours = [true false];
    labelCones = false;
    subplotTitle = '';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 2;
    iteration = 1;
    labelContours = [false false];
    labelCones = false;
    subplotTitle = sprintf('%04.0f', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, false, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 3;
    iteration = intermediateIterationsToDisplay(1);
    labelContours = [false false];
    labelCones = false;
    subplotTitle =  sprintf('%04.0f', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 4;
    iteration = intermediateIterationsToDisplay(2);
    labelContours = [false false];
    labelCones = false;
    subplotTitle =  sprintf('%04.0f', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 5;
    iteration = size(obj.latticeAdjustmentSteps, 1);
    labelContours = [false false];
    labelCones = false;
    subplotTitle =  sprintf('%04.0f', iteration);
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 6;
    iteration = size(obj.latticeAdjustmentSteps, 1);
    labelContours = [false false];
    labelCones = true;
    contourLevels = [];
    subplotTitle = '';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    % Restore the final state of coneLocsHexGrid
    obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(end, :, :));
end

function plotMosaic(obj, subplotPosVectors, row, col, ...
    displayedXrangeMicrons, displayedYrangeMicrons, plotHexMesh, ...
    plotContours, contourLevels, labelContours, labelCones, ...
    backgroundColor, iteration, subplotTitle)
% Plot the mosaic
%
% Syntax:
%   plotMosaic(obj, subplotPosVectors, row, col, ...
%       displayedXrangeMicrons, displayedYramgeMicrons, plotHexMesh, ...
%       plotContours, contourLevels, labelContours, labelCones, ...
%       backgroundColor, iteration, subplotTitle)
%
% Description:
%    Plot the mosaic
%
% Inputs:
%    obj                    - The cone mosaic hex object
%    subplotPosVectors      - A matrix struct containing vectors for
%                             subplot positions.
%    row                    - The position row
%    col                    - The position column
%    displayedXrangeMicrons - The X range to display, in microns
%    displayedYramgeMicrons - The Y range to display, in microns
%    plotHexMesh            - The plot for the hex mesh
%    plotContours           - The plot for the contours
%    contourLevels          - The number of contour levels
%    labelContours          - The contour labels
%    labelCones             - The cone labels
%    backgroundColor        - The color to use for background
%    iteration              - The iteration for cone locations
%    subplotTitle           - The title for the subplot
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    % Get cone locs at desired iteration
    if (iteration == 0)
        obj.coneLocsHexGrid = obj.initialLattice;
    else
        obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(...
            iteration, :, :));
    end

        
    ax = axes('Position', subplotPosVectors(row, col).v, ...
        'units', 'normalized', 'Color', backgroundColor);


    if (isempty(displayedXrangeMicrons))
        displayedXrangeMeters = obj.center(1)+ obj.width / 2 * [-1 1];
    else
        displayedXrangeMeters = displayedXrangeMicrons / 2 * 1e-6 * [-1 1];
    end
    if (isempty(displayedYrangeMicrons))
        displayedYrangeMeters = obj.center(2)+ [0 obj.height / 2];
    else
        displayedYrangeMeters = [0 displayedYrangeMicrons / 2 * 1e-6];
    end

    % Determine cones to be plotted
    coneApertureMicrons = obj.lambdaMin / 1e6;
    kFactorX = -1.5;
    kFactorY = -1.5;

    idx = find(...
            (obj.coneLocsHexGrid(:, 1) - kFactorX * coneApertureMicrons ...
            > displayedXrangeMeters(1)) & ...
            (obj.coneLocsHexGrid(:, 1) + kFactorX * coneApertureMicrons ...
            < displayedXrangeMeters(2)) & ...
            (obj.coneLocsHexGrid(:, 2) - kFactorY * coneApertureMicrons ...
            > displayedYrangeMeters(1)) & ...
            (obj.coneLocsHexGrid(:, 2) + kFactorY * coneApertureMicrons ...
            < displayedYrangeMeters(2)));

    if (labelCones)
        obj.visualizeGrid('axesHandle', ax, ...
            'visualizedConeAperture', 'both', ...
            'backgroundColor', backgroundColor);
    else
        % Plot the cones
        edgeColor = [0.1 0.1 0.1];
        faceColor = 0.99*[1 1 1];
        lineStyle = '-';
        lineWidth = 0.5;
        iTheta = (0:30:360) / 180 * pi;
        coneApertureRadius = obj.lambdaMin/2;
        coneAperture.x = coneApertureRadius*cos(iTheta) * 1e-6;
        coneAperture.y = coneApertureRadius*sin(iTheta) * 1e-6;
        coneMosaicHex.renderPatchArray(ax, coneAperture, ...
            squeeze(obj.coneLocsHexGrid(idx,1)), ...
            squeeze(obj.coneLocsHexGrid(idx,2)), ...
            edgeColor, faceColor, lineStyle, lineWidth);

    end
    hold on;

    % superimpose density contours
    if (plotContours > 0) && (~isempty(contourLevels))
        % Compute actual mosaic density
        [densityMapMosaic, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from mosaic');

        
        % Compute model mosaic density
        [densityMapModel, densityMapSupportX, densityMapSupportY] = ...
            obj.computeDensityMap('from model');

        
        contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapModel, contourLevels, 'LineColor', [0.7 0.7 1.0], ...
            'LineWidth', 4.0, 'ShowText', 'off');
        
        
        [cH1, hH1] = contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapModel, contourLevels, 'LineColor', [0.0 .0 1.0], ...
            'LineWidth', 2, 'ShowText', 'off');
        
        contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapMosaic, contourLevels, 'LineColor', [1 0.7 0.7], ...
            'LineWidth', 4.0, 'ShowText', 'off');
        
        [cH2, hH2] = contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapMosaic, contourLevels, 'LineColor', [1 0.0 0.0], ...
            'LineWidth', 2, 'ShowText', 'off');
        
        if (labelContours(1))
            clabel(cH1, hH1, 'FontWeight', 'bold', 'FontName', 'Menlo', ...
                'FontSize', 14, 'Color', [0.1 0.1 1.0], ...
                'Background', [1 1 1], 'LabelSpacing', 400);
        end
        
        if (labelContours(2))
            clabel(cH2, hH2, 'FontWeight', 'bold', 'FontName', 'Menlo', ...
                'FontSize', 14, 'Color', [1 0.1 0.1], ...
                'Background', [1 1 1], 'LabelSpacing', 200);
        end
    end

    if (plotHexMesh)
        coneMosaicHex.renderHexMesh(ax, obj.coneLocsHexGrid(idx, 1), ...
            obj.coneLocsHexGrid(idx, 2), [0.5 0.5 0.5], 'none', ...
            0.8, 0.2, '-');
    end

    xRange = displayedXrangeMeters + [-2 2]*1e-6;
    yRange = displayedYrangeMeters + [-2 2]*1e-6;
    
    % Plot box around plot.
    xBox = [xRange(1) xRange(1) xRange(2) xRange(2) xRange(1)];
    yBox = [yRange(1) yRange(2) yRange(2) yRange(1) yRange(1)];
    plot(ax, xBox, yBox, 'k-', 'LineWidth', 1.0);
    
    
    % Finalize plot
    hold off;
    axis 'equal';
    axis 'xy';
    xTicks = 1e-6 * (-1000:20:1000);
    xTickLabels = sprintf('%2.0f\n', xTicks / 1e-6);
    yTicks = 1e-6 * (-1000:20:1000);
    yTickLabels = sprintf('%2.0f\n', yTicks / 1e-6);

    yAxisColor = 'k';
    xAxisColor = 'k';


    set(ax, 'XLim', xRange, 'YLim', yRange, ...
        'XTick', xTicks, 'XTickLabels', xTickLabels, ...
        'YTick', yTicks, 'YTickLabels', yTickLabels, ...
        'FontSize', 18, 'XColor', xAxisColor, 'YColor', yAxisColor, ...
        'LineWidth', 1.0);
    grid off;
    box off;
    
    if (row == size(subplotPosVectors, 1))
        xlabel('\it microns', 'Interpreter','tex', 'FontWeight', 'normal', 'FontSize', 24);
    else
        xAxisColor = 'none';
        set(ax, 'xTickLabels', {});
    end

    if (col == 1)
        ylabel('\it microns', 'Interpreter','tex', 'FontWeight', 'normal', 'FontSize', 24); 
    else
        yAxisColor = 'none'; 
    end
    
    if (~isempty(subplotTitle))
        t = text(64*1e-6, 61*1e-6 , subplotTitle, 'Color', [0 0 0], 'FontSize', 16);
        t.BackgroundColor = 0.9*[1 1 1];
        t.EdgeColor = 0.3*[1 1 1];
    end
    drawnow;
end