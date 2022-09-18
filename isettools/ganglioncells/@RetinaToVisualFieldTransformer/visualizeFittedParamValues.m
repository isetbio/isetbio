function visualizeFittedParamValues(params)
    n = numel(params.initialValues);
    if (n < 5)
        rows = 1;
    else
        rows = ceil(n/5);
    end

    hFig = figure(999); clf;
    set(hFig, 'Position', [10 10 1600 500], 'Color', [1 1 1 ]);
    for iParam = 1:n
        ax = subplot(rows,5,iParam);
        visualizeParamRange(ax, ...
                params.names{iParam}, ...
                params.scaling{iParam}, ...
                params.initialValues(iParam), ...
                params.finalValues(iParam), ...
                params.lowerBounds(iParam), ...
                params.upperBounds(iParam) ...
                );
    end
end

function visualizeParamRange(ax, paramName, paramScaling, initialValue, finalValue, lowerBound, upperBound)
    
    plot(ax, [0 0], [lowerBound upperBound], 'k-', 'LineWidth', 6, 'Color', [0.8 0.8 0.7]); hold(ax, 'on');
    plot(ax, [0 0], [initialValue finalValue], 'k--', 'LineWidth', 1.0);
    plot(ax, 0, initialValue, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', [1 0.5 0.50]);
    plot(ax, 0, finalValue, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.50], 'LineWidth', 1.5);
    axis(ax, 'square');
    if (strcmp(paramScaling, 'log'))
        ticks = logspace(log10(lowerBound), log10(upperBound), 5);
    else
        ticks = linspace(lowerBound, upperBound, 5);
    end
    grid(ax, 'on');
    if (lowerBound == upperBound)
        ticks = lowerBound;
        if (strcmp(paramScaling, 'linear'))
            lowerBound = lowerBound-0.5;
            upperBound = upperBound+0.5;
        else
            lowerBound = lowerBound/2;
            upperBound = upperBound*2;
        end
    end

    set(ax, 'XLim', [-1 1], 'XTick', [], 'YTick', ticks, 'YLim', [lowerBound upperBound], ...
        'YScale', paramScaling);
    set(ax, 'FontSize', 16);
    xlabel(paramName);
end

