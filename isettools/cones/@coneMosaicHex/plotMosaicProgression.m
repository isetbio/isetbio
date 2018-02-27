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
%    **Needed!**
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
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 880 1300]);

    rowsNum = 6;
    colsNum = 1;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.015, ...
           'widthMargin', 0.00, ...
           'leftMargin', 0.02, ...
           'rightMargin', -0.02, ...
           'bottomMargin', 0.05, ...
           'topMargin', 0.00);

    backgroundColor = [1 1 1];

    col = 1;
    row = 1;
    iteration = 0;
    labelContours = false;
    labelCones = false;
    subplotTitle = '(A)';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, false, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 2;
    iteration = 1;
    labelContours = false;
    labelCones = false;
    subplotTitle = '(B)';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, false, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 3;
    iteration = intermediateIterationsToDisplay(1);
    labelContours = false;
    labelCones = false;
    subplotTitle = '(C)';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 4;
    iteration = intermediateIterationsToDisplay(2);
    labelContours = false;
    labelCones = false;
    subplotTitle = '(D)';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 5;
    iteration = size(obj.latticeAdjustmentSteps, 1);
    labelContours = true;
    labelCones = false;
    subplotTitle = '(E)';
    plotMosaic(obj, subplotPosVectors, row, col, ...
        displayedXrangeMicrons, displayedYrangeMicrons, false, true, ...
        contourLevels, labelContours, labelCones, backgroundColor, ...
        iteration, subplotTitle);

    row = 6;
    iteration = size(obj.latticeAdjustmentSteps, 1);
    labelContours = false;
    labelCones = true;
    contourLevels = [];
    subplotTitle = '(F)';
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
    kFactorX = 0.5;
    kFactorY = -0.5;
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
            'visualizedConeAperture', 'both');
    else
        % Plot the cones
        edgeColor = 'none';
        faceColor = 0.7*[1 1 1];
        lineStyle = '-';
        iTheta = (0:60:300) / 180 * pi;
        coneApertureRadius = obj.lambdaMin / 2;
        coneAperture.x = coneApertureRadius * cos(iTheta) * 1e-6;
        coneAperture.y = coneApertureRadius * sin(iTheta) * 1e-6;
        coneMosaicHex.renderPatchArray(ax, coneAperture, ...
            squeeze(obj.coneLocsHexGrid(idx, 1)), ...
            squeeze(obj.coneLocsHexGrid(idx, 2)), ...
            edgeColor, faceColor, lineStyle);
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
   
        [cH, hH] = contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapModel, contourLevels, 'LineColor', [0.1 .1 1.0], ...
            'LineWidth', 3.0, 'ShowText', 'off');
        if (labelContours)
            clabel(cH, hH, 'FontWeight', 'bold', 'FontName', 'Menlo', ...
                'FontSize', 12, 'Color', [0.1 0.1 1.0], ...
                'Background', [1 1 1], 'LabelSpacing', 400);
        end

        [cH, hH] = contour(ax, densityMapSupportX, densityMapSupportY, ...
            densityMapMosaic, contourLevels, 'LineColor', [1 0.1 0.1], ...
            'LineWidth', 3.0, 'ShowText', 'off');
        if (labelContours)
            clabel(cH, hH, 'FontWeight', 'bold', 'FontName', 'Menlo', ...
                'FontSize', 12, 'Color', [1 0.1 0.1], ...
                'Background', [1 1 1], 'LabelSpacing', 200);
        end
    end

    if (plotHexMesh)
        coneMosaicHex.renderHexMesh(ax, obj.coneLocsHexGrid(idx, 1), ...
            obj.coneLocsHexGrid(idx, 2), [0.5 0.5 0.5], 'none', ...
            0.8, 0.2, '-');
    end

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

    if (row == size(subplotPosVectors, 1))
        xlabel('microns');
    else
        xAxisColor = 'none';
    end

    if (col == 1), ylabel('microns'); else, yAxisColor = 'none'; end
    
    set(gca, ...
        'XLim', [displayedXrangeMeters(1) - 2 * 1e-6, ...
        displayedXrangeMeters(2) + 3 * 1e-6], ...
        'YLim', [displayedYrangeMeters(1) - 2 * 1e-6, ...
        displayedYrangeMeters(2) + 3 * 1e-6], ...
        'XTick', xTicks, 'XTickLabels', xTickLabels, ...
        'YTick', yTicks, 'YTickLabels', yTickLabels, ...
        'TickDir', 'both', ...
        'FontSize', 18, ...
        'XColor', xAxisColor, ...
        'YColor', yAxisColor, ...
        'LineWidth', 1.0);
    grid off;
    box off;

    text(displayedXrangeMeters(1) - 40 * 1e-6, ...
        displayedYrangeMeters(2) - 2 * 1e-6, subplotTitle, ...
        'Color', [0 0 0], 'FontSize', 26, 'FontWeight', 'bold');
    drawnow;
end