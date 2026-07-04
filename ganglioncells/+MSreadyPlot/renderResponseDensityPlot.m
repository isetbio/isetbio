function renderResponseDensityPlot(ax, x, y, response, ...
    responseRange, responseColorMap, ...
    xLabelTitle, yLabelTitle, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    % Plot
    imagesc(ax,x,y, response);

    % xy axis labels
    if (~noXLabel)
        xlabel(ax, xLabelTitle, 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, yLabelTitle, 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', {});
    end

    % ticks and grids
    grid(ax, 'on')
    xtickangle(ax, 0);
    set(ax, 'XLim', [x(1) x(end)], 'YLim', [y(1) y(end)], 'CLim', responseRange);
    set(ax, 'XTick', [], 'YTick', []);
    set(ax, 'TickDir', 'both');

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % Colormap
    colormap(responseColorMap);

    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end