function renderSTF(ax, sfSupportCPD, compositeSTF, compositeSTFfit, centerSTFfit, surroundSTFfit, plotTitle, theLegends, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    
    maxAll = max([max(compositeSTF(:)) max(compositeSTFfit(:)) max(centerSTFfit(:)) max(surroundSTFfit(:))]);
    compositeSTF = compositeSTF / maxAll;
    compositeSTFfit = compositeSTFfit / maxAll;
    centerSTFfit = centerSTFfit / maxAll;
    surroundSTFfit = surroundSTFfit / maxAll;

    y1 = centerSTFfit;
    y2 = 0*centerSTFfit;
    faceColor = [1 0.5 0.5];
    edgeColor = [1 0 0 ];

    faceAlpha = 0.3;
    lineWidth = 1.0;
    lineStyle = '-';
    p3 = shadedAreaBetweenTwoLines(ax, sfSupportCPD, y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
    
    hold(ax, 'on');

    y1 = surroundSTFfit;
    y2 = 0*surroundSTFfit;
    faceColor = [0.5 0.5 1];
    edgeColor = [0 0 1];

    faceAlpha = 0.3;
    lineWidth = 1.0;
    lineStyle = '-';
    p4 = shadedAreaBetweenTwoLines(ax, sfSupportCPD, y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);


    p2 = plot(ax, sfSupportCPD, compositeSTFfit, 'k:', 'LineWidth', ff.lineWidth*2);
    p1 = plot(ax, sfSupportCPD, compositeSTF, 'go', ...
        'MarkerEdgeColor', [0.3 0.5 0.3], 'MarkerFaceColor', [0.6 0.8 0.8], 'MarkerSize', ff.markerSize-2, ...
        'LineWidth', ff.lineWidth*1.5);
    axis(ax, 'square');
    
    % ticks and grids
    xTicks = [0.1 0.3 1 3 10 30 100];
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', [min(sfSupportCPD)-0.05 100], 'YLim', [-0.02 1]);
    set(ax, 'XTick', xTicks, 'YTick', 0:0.2:1, ...
            'XTickLabel', xTicks, 'YTickLabel', 0:0.2:1);
    set(ax, 'TickDir', 'both')
    set(ax, 'XScale', 'log');

    % xy axis labels
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    if (~noXLabel)
        xlabel(ax, 'spatial frequency (c/deg)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'STF', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', 0:0.2:1, 'YTickLabel', sprintf('%1.1f\n', 0:0.2:1));
    end

    % legends
    if (~isempty(theLegends))
        hl = legend(ax, [p1 p2 p3 p4], theLegends, ...
            'Location', 'SouthWest', 'NumColumns', 1, ...
            'FontSize', ff.legendFontSize, 'Box', 'off');
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

 function p = shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
 end