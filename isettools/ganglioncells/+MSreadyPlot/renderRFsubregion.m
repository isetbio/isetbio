function renderRFsubregion(ax, xSupportDegs, ySupportDegs, theSubregionMap, rfRangeDegs, plotTitle, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('withLineWeightingFunction', [], @isnumeric);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    lineWeightingFunction = p.Results.withLineWeightingFunction;
    
    cmap = brewermap(1024, 'greys');
    
    contourPlot = false;
    if (contourPlot)
        zLevels = [0.1 0.3 0.5 0.7 0.9];
        MSreadyPlot.transparentContourPlot(ax, ...
            xSupportDegs*60, ySupportDegs*60, ...
            theSubregionMap/max(theSubregionMap(:)), ...
            zLevels, cmap, ...
            'faceAlpha', 0.8);
    else
        imagesc(ax, xSupportDegs*60, ySupportDegs*60, theSubregionMap/max(theSubregionMap(:)));
        colormap(ax, cmap)
    end

    if (~isempty(lineWeightingFunction))
        hold(ax, 'on');
        rfRangeArcMin = rfRangeDegs*60;
        baseline = -rfRangeArcMin;
        amplitude = 1.9*rfRangeArcMin;
        
        y2 = baseline + 0*lineWeightingFunction;
        if (max(lineWeightingFunction) < abs(min(lineWeightingFunction)))
            faceColor = [0.5 0.5 1];
            edgeColor = [0 0 1];
            y1 = baseline - amplitude*lineWeightingFunction;
        else
            y1 = baseline + amplitude*lineWeightingFunction;
            faceColor = [1 0.5 0.5];
            edgeColor = [1 0 0 ];
        end
        xArcMin = 60*xSupportDegs;

        % Only show down to 1% of the peak
        thresholdNormalizedValue = 1/100;
        
        idx = find(abs(lineWeightingFunction)>= thresholdNormalizedValue * max(abs(lineWeightingFunction)));
        xArcMin = xArcMin(idx);
        y1 = y1(idx);
        y2 = y2(idx);

        faceAlpha = 0.3;
        lineWidth = 1.0;
        lineStyle = '-';
        shadedAreaBetweenTwoLines(ax, xArcMin', y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
    end

    set(ax, 'CLim', [0 1], 'ZLim', [0 1]);
    axis(ax, 'image');
    axis(ax, 'xy');
    

    rfRangeArcMin = rfRangeDegs*60;
    rfTicksArcMin = MSreadyPlot.spatialMapTicksArcMin(rfRangeArcMin);

    if (rfRangeArcMin <= 2)
        rfTickLabels = sprintf('%2.1f\n', rfTicksArcMin);
    else
        rfTickLabels = sprintf('%2.0f\n', rfTicksArcMin);
    end

    % ticks and grids
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'XLim', rfRangeArcMin*1.05*[-1 1], 'YLim', rfRangeArcMin*1.05*[-1 1]);
    set(ax, 'XTick', rfTicksArcMin , 'YTick', rfTicksArcMin, ...
            'XTickLabel', rfTickLabels, 'YTickLabel', rfTickLabels);
    set(ax, 'TickDir', 'in')
    
    % xy axis labels
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    if (~noXLabel)
        xlabel(ax, 'space, x (arc min)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'space, y (arc min)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', {});
    end

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end

 function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
 end
