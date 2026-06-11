function renderSTF(ax, sfSupportCPD, compositeSTF, sfSupportCDPfit, compositeSTFfit, centerSTFfit, surroundSTFfit, plotTitle, theLegends, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTickLabel', false, @islogical);
    p.addParameter('noYTickLabel', false, @islogical);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('visualizeCronerKaplanTypicalSTF', false, @(x)(islogical(x)||(isnumeric(x))));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTickLabel = p.Results.noXTickLabel;
    noYTickLabel = p.Results.noYTickLabel;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    visualizeCronerKaplanTypicalSTF = p.Results.visualizeCronerKaplanTypicalSTF;

    maxAll = max([max(compositeSTF(:)) max(compositeSTFfit(:)) max(centerSTFfit(:)) max(surroundSTFfit(:))]);
    compositeSTF = compositeSTF / maxAll;
    compositeSTFfit = compositeSTFfit / maxAll;
    centerSTFfit = centerSTFfit / maxAll;
    surroundSTFfit = surroundSTFfit / maxAll;


    minSF = min(sfSupportCDPfit);
    if (~isempty(visualizedSpatialFrequencyRange))
        minSF = visualizedSpatialFrequencyRange(1);
    end

    idx = find(sfSupportCDPfit >= minSF);
    sfSupportCDPfit = sfSupportCDPfit(idx);
    compositeSTFfit = compositeSTFfit(idx);
    surroundSTFfit = surroundSTFfit(idx);
    centerSTFfit = centerSTFfit(idx);

    y1 = centerSTFfit;
    y2 = 0*centerSTFfit;
    faceColor = [0.75 0.75 0.75];
    edgeColor = [0.3 0.3 0.3];

    faceAlpha = 0.3;
    lineWidth = 1.0;
    lineStyle = '-';
    p3 = shadedAreaBetweenTwoLines(ax, sfSupportCDPfit, y1, y2, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
    
    hold(ax, 'on');

    y1 = surroundSTFfit;
    y2 = 0*surroundSTFfit;
    faceColor = [0.4 0.4 0.4];
    edgeColor = [0.2 0.2 0.2];

    faceAlpha = 0.3;
    lineWidth = 1.0;
    lineStyle = '-';
    p4 = shadedAreaBetweenTwoLines(ax, sfSupportCDPfit, y1, y2, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);


    
    p2 = plot(ax, sfSupportCDPfit, compositeSTFfit, 'r-', 'LineWidth', ff.lineWidth*2);
    p1 = plot(ax, sfSupportCPD, compositeSTF, 'ro', ...
        'MarkerEdgeColor', [1.0 0 0], 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', ff.markerSize, ...
        'LineWidth', ff.lineWidth);

    if (islogical(visualizeCronerKaplanTypicalSTF))
        if (visualizeCronerKaplanTypicalSTF)
        typicalSTF = RGCmodels.CronerKaplan.digitizedData.typicalRGCSTF();
        plot(ax, typicalSTF.sfCPD, typicalSTF.amplitude*max(compositeSTFfit), 'b--', 'LineWidth', ff.lineWidth);
        end
    else
        typicalSTF = visualizeCronerKaplanTypicalSTF;
        plot(ax, sfSupportCPD, typicalSTF/max(typicalSTF)*max(compositeSTFfit), 'b--', 'LineWidth', ff.lineWidth);
    end

    axis(ax, 'square');
    
    % ticks and grids
    xTicks = [0.01 0.03 0.1 0.3 1 3 10 30 100];
    if (isempty(visualizedSpatialFrequencyRange))
         visualizedSpatialFrequencyRange(1) = minSF;
        visualizedSpatialFrequencyRange(2) = 100;
    end
    
    visualizedSpatialFrequencyRange = [visualizedSpatialFrequencyRange(1) visualizedSpatialFrequencyRange(2)];
    
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', visualizedSpatialFrequencyRange, 'YLim', [1*ff.axisOffsetFactor 1]);
    
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

    if (~noYTickLabel)
    else
        set(ax, 'YTickLabel', {});
    end

    if (~noXTickLabel)
    else
        set(ax, 'XTickLabel', {});
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

 function p = shadedAreaBetweenTwoLines(ax,x,y1, y2, ...
     faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
 end