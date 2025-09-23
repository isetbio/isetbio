%
% RGCMosaicAnalyzer.visualize.setLimsAndTicks()
%

function setLimsAndTicks(ax, domainVisualizationTicks, domainVisualizationLimits)

    minTickIncrement = round(min([min(abs(diff(domainVisualizationTicks.x))) min(abs(diff(domainVisualizationTicks.y)))])*10)/10;

    if (minTickIncrement >= 1)
        set(ax, 'XTickLabel', sprintf('%1.0f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.0f\n', domainVisualizationTicks.y));
    elseif (minTickIncrement >= 0.1)
        set(ax, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
    elseif (minTickIncrement >= 0.01)
        set(ax, 'XTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.y));
    else
        set(ax, 'XTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.y));
    end
    set(ax, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));
    set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);

    %xlabel(ax, 'eccentricity, x (degs)');
    %ylabel(ax, 'eccentricity, y (degs)');
end