function renderConeFundamentals(ax, wavelengthSupport, Lcone, Mcone, Scone, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    % Plot
    hold(ax, 'on');
    legends = {};
    if (~isempty(Lcone))
        scatter(ax, wavelengthSupport, Lcone, (ff.markerSize)^2, ...
            'MarkerFaceColor', ff.LconeColor,  'MarkerEdgeColor', ff.LconeColor*0.5,  ...
            'MarkerFaceAlpha', 0.6, 'LineWidth', ff.lineWidth);
        legends{numel(legends)+1} = 'L';
    end

    if (~isempty(Mcone))
        scatter(ax, wavelengthSupport, Mcone, (ff.markerSize)^2, ...
            'MarkerFaceColor', ff.MconeColor,  'MarkerEdgeColor', ff.MconeColor*0.5,  ...
            'MarkerFaceAlpha', 0.6, 'LineWidth', ff.lineWidth);
        legends{numel(legends)+1} = 'M';
    end

    if (~isempty(Scone))
        scatter(ax, wavelengthSupport, Scone, (ff.markerSize)^2, ...
            'MarkerFaceColor', ff.SconeColor,  'MarkerEdgeColor', ff.SconeColor*0.5,  ...
            'MarkerFaceAlpha', 0.6, 'LineWidth', ff.lineWidth);
        legends{numel(legends)+1} = 'S';
    end

    % axis(ax, 'square');

    % xy axis labels
    if (~noXLabel)
        xlabel(ax, 'wavelength', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'weight', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', {});
    end

    % ticks and grids
    grid(ax, 'on')
    xtickangle(ax, 0);
    set(ax, 'XLim', [wavelengthSupport(1)-10 wavelengthSupport(end)], 'YLim', [-0.02 0.5]);
    set(ax, 'XTick', 400:50:700, 'YTick', 0:0.1:1);
    set(ax, 'TickDir', 'both');

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % Legends
    legend(ax, legends)

    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end