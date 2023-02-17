function renderRFcompositeMap(ax, xSupportDegs, ySupportDegs, theCompositeRFmap, rfRangeDegs, plotTitle, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    
    
    compositeMapLineWeightingFunctionX = sum(theCompositeRFmap,1);
    compositeMapLineWeightingFunctionX = compositeMapLineWeightingFunctionX / max(abs(compositeMapLineWeightingFunctionX));

    sensitivityRange = [-0.35 1.05];

    plot(ax, xSupportDegs*60, compositeMapLineWeightingFunctionX, 'k-', 'LineWidth', 1.5);
    
    
    set(ax, 'CLim', [0 1]);
    axis(ax, 'square');
    

    rfRangeArcMin = rfRangeDegs*60;
    rfTicksArcMin = MSreadyPlot.spatialMapTicksArcMin(rfRangeArcMin);

    if (rfRangeArcMin <= 2)
        rfTickLabels = sprintf('%2.1f\n', rfTicksArcMin);
    else
        rfTickLabels = sprintf('%2.0f\n', rfTicksArcMin);
    end

    % ticks and grids
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'XLim', rfRangeArcMin*1.05*[-1 1], 'YLim', sensitivityRange);
    set(ax, 'XTick', rfTicksArcMin , 'YTick', -1:0.2:1, ...
            'XTickLabel', rfTickLabels, 'YTickLabel', -1:0.2:1);
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
