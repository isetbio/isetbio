function renderResponseTimeSeriesPlot(ax, x, response, ...
    responseRange, responseColorMap, ...
    xLabelTitle, yLabelTitle, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('MarkerEdgeColor', [0 0 0], @(x)(isnumeric(x)&&(numel(x)==3)));
    p.addParameter('MarkerFaceColor', [0.7 0.7 0.7], @(x)(isnumeric(x)&&(numel(x)==3)));
    p.addParameter('markerSize', 10, @isscalar);
    p.addParameter('lineWidth', 1.5, @isscalar);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    MarkerFaceColor = p.Results.MarkerFaceColor;
    MarkerEdgeColor = p.Results.MarkerEdgeColor;
    markerSize = p.Results.markerSize;
    lineWidth = p.Results.lineWidth;

    % Plot
    hold(ax, 'on');
    responseInstances = size(response,1);

    if (responseInstances == 1)
        iInstance = 1;
        plot(ax,x, response(iInstance,:), 'ko-', 'LineWidth', lineWidth, 'Color', MarkerEdgeColor, ...
                'MarkerSize', markerSize, 'MarkerEdgeColor', MarkerEdgeColor, 'MarkerFaceColor', MarkerFaceColor);
    else
        for iInstance = 1:responseInstances
            plot(ax,x, response(iInstance,:), 'k.-', 'LineWidth', lineWidth, 'Color', MarkerEdgeColor)
        end
    end

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
    yTicks = [responseRange(1) 0 responseRange(2)];
    set(ax, 'XLim', [x(1) x(end)], 'YLim', responseRange);
    set(ax, 'XTick', x, 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
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