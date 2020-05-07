function plotSensitivities(theAxes, d, model, pointSize, color,displayYLabel, theLabel)

    dSens = d('sensitivity');
    dSize = d('size');
    
    % Face color
    c = color + 0.5*[1 1 1];
    c(c>1)=1;
    
    % Plot the digitized data
    hold(theAxes, 'on');
    scatter(theAxes, dSens.radiusDegs, dSens.peakSensitivity, pointSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', color);
    if (isfield(dSize, 'radiusDegsMedianTable'))
        scatter(theAxes, dSize.radiusDegsMedianTable, dSens.peakSensitivityMedianTable, 's', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', [0.2 0.2 0.2]);
    end
    
     % Plot the fitted model
    plot(theAxes, model.radiusDegs, model.function(model.params, model.radiusDegs), 'k-');
    
    axis(theAxes, 'square');
    set(theAxes, 'XLim', [0.01 3], 'YLim', [0.003 3000]);
    set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 1 3], 'XTickLabel', {0.01 0.03 0.1 0.3 1 3}, ...
                 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000], ...
                 'YTickLabel', {'0.003' '0.01' '0.03' '0.1' '0.3' '1' '3' '10' '30' '100' '300' '1000', '3000'});
    xlabel(theAxes, 'radius (degs)');
    if (displayYLabel)
        ylabel(theAxes, 'peak sensitivity (imp/sec/contrast/deg^2)');
    end
    if (~isempty(theLabel))
    	title(theAxes, theLabel, 'Color', color, 'FontSize', 30);
    end
    set(theAxes, 'XScale', 'log', 'YScale', 'log');
    grid(theAxes, 'on');
    box(theAxes, 'off');
end