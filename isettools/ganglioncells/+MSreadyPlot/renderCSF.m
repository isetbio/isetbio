function renderCSF(ax, sfSupportCPD, CSF, plotTitle, theLegends, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTickLabel', false, @islogical);
    p.addParameter('noYTickLabel', false, @islogical);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('visualizedSensitivityRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('colors', []);
    p.addParameter('markers', []);

    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTickLabel = p.Results.noXTickLabel;
    noYTickLabel = p.Results.noYTickLabel;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    visualizedSensitivityRange = p.Results.visualizedSensitivityRange;
    colors = p.Results.colors;
    markers = p.Results.markers;
    if (isempty(colors))
        colors = [1 0 0; 0 0 1; 0 0 0; 0 0.8 0];
    end
    if (isempty(markers))
        for iPlot = 1:size(CSF,1)
            markers{iPlot} = 'o';
        end
    end

    plotHanlesList = [];
    

    for iPlot = 1:size(CSF,1)
        theColor = colors(iPlot,:);
        plotHanlesList(iPlot) = plot(ax, sfSupportCPD, CSF(iPlot,:), 'o-', ...
            'Marker', markers{iPlot}, 'MarkerEdgeColor', theColor, 'MarkerFaceColor', 0.5*theColor + 0.5*[1 1 1], ...
            'MarkerSize', ff.markerSize*1.25, ...
            'Color', theColor, ...
            'LineWidth', ff.lineWidth*1.5);
        hold(ax, 'on');
    end

    % ticks and grids
    xTicks = [0.01 0.03 0.1 0.3 1 3 10 30 100];

    if (isempty(visualizedSpatialFrequencyRange))
        visualizedSpatialFrequencyRange(1) = 0.01;
        visualizedSpatialFrequencyRange(2) = 100;
    end
    
    if (isempty(visualizedSensitivityRange))
        visualizedSensitivityRange = [0 max(CSF)];
    end

    visualizedSpatialFrequencyRange = [visualizedSpatialFrequencyRange(1) visualizedSpatialFrequencyRange(2)];
    
    grid(ax, 'on'); box(ax, 'off');
    set(ax, ...
        'XLim', [visualizedSpatialFrequencyRange(1) visualizedSpatialFrequencyRange(2)], ...
        'YLim', [visualizedSensitivityRange(1)+(visualizedSensitivityRange(2)-visualizedSensitivityRange(1))*ff.axisOffsetFactor visualizedSensitivityRange(2)]);
    
    set(ax, 'XTick', xTicks, 'YTick', 0:100:1000, ...
            'XTickLabel', sprintf('%.1f\n',xTicks));
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
        ylabel(ax, 'contrast sensitivity', 'FontAngle', ff.axisFontAngle);
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
        hl = legend(ax, plotHanlesList, theLegends, ...
            'Location', 'SouthWest', 'NumColumns', 1, ...
            'FontSize', ff.legendFontSize, 'Box', 'off', 'Color', [1 1 1]);
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