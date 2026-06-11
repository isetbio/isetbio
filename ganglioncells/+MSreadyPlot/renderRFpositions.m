function renderRFpositions(ax, xyCoords, xLimsDegs, yLimsDegs, plotTitle, ff, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});

    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    plot(ax, xyCoords(:,1), xyCoords(:,2), ...
         'k.', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', ff.markerSize/2);

    axis(ax, 'equal');
    set(ax, 'XLim', xLimsDegs, 'YLim', yLimsDegs);

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