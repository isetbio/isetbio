function renderCSF(ax, sfSupportCPD, CSF, plotTitle, theLegends, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTickLabel', false, @islogical);
    p.addParameter('noYTickLabel', false, @islogical);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTickLabel = p.Results.noXTickLabel;
    noYTickLabel = p.Results.noYTickLabel;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;

    p2 = plot(ax, sfSupportCPD, CSF, 'ro-', ...
            'MarkerEdgeColor', [1.0 0 0], 'MarkerFaceColor', [1 0.5 0.5], ...
            'MarkerSize', ff.markerSize, ...
            'LineWidth', ff.lineWidth);
    axis(ax, 'square');
    
    % ticks and grids
    xTicks = [0.01 0.03 0.1 0.3 1 3 10 30 100];
    if (isempty(visualizedSpatialFrequencyRange))
        visualizedSpatialFrequencyRange(1) = minSF;
        visualizedSpatialFrequencyRange(2) = 100;
    end
    
    visualizedSpatialFrequencyRange = [visualizedSpatialFrequencyRange(1) visualizedSpatialFrequencyRange(2)];
    
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', visualizedSpatialFrequencyRange, 'YLim', [100*ff.axisOffsetFactor 100]);
    
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
        ylabel(ax, 'sensitivity', 'FontAngle', ff.axisFontAngle);
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