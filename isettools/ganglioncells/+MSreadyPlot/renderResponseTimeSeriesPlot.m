function renderResponseTimeSeriesPlot(ax, temporalSupport, response, ...
    responseRange, ...
    xLabelTitle, yLabelTitle, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('xTicks', [], @isnumeric);
    p.addParameter('yTicks', [], @isnumeric);
    p.addParameter('lineColor', [0 0 0], @(x)(isnumeric(x)&&(numel(x)==3)));
    p.addParameter('markerEdgeColor', [0 0 0], @(x)(isnumeric(x)&&(numel(x)==3)));
    p.addParameter('markerFaceColor', [0.7 0.7 0.7], @(x)(isnumeric(x)&&(numel(x)==3)));
    p.addParameter('markerSize', 10, @isscalar);
    p.addParameter('lineWidth', 1.5, @isscalar);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    lineColor = p.Results.lineColor;
    markerFaceColor = p.Results.markerFaceColor;
    markerEdgeColor = p.Results.markerEdgeColor;
    markerSize = p.Results.markerSize;
    lineWidth = p.Results.lineWidth;
    xTicks = p.Results.xTicks;
    yTicks = p.Results.yTicks;

    % Plot
    hold(ax, 'on');
    responseInstances = size(response,1);

    if (responseInstances == 1)
        iInstance = 1;
        plot(ax, temporalSupport, response(iInstance,:), 'ko-', 'LineWidth', lineWidth, 'Color', lineColor, ...
                'MarkerSize', markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor);
    else
        for iInstance = 1:responseInstances
            plot(ax, temporalSupport, response(iInstance,:), '-', 'LineWidth', lineWidth, 'Color', lineColor)
        end
    end

    

    % ticks and grids
    grid(ax, 'on')
    xtickangle(ax, 0);
    
    dt = ff.axisOffsetFactor * (temporalSupport(end)-temporalSupport(1));
    dR = ff.axisOffsetFactor * (responseRange(end)-responseRange(1));
    set(ax, 'XLim', [temporalSupport(1)+dt  temporalSupport(end)], 'YLim', [responseRange(1)+dR responseRange(2)]);
    if (~isempty(xTicks))
        set(ax, 'XTick', xTicks);
    end
    if (~isempty(yTicks))
        set(ax, 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
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
    
    set(ax, 'TickDir', 'both');

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end