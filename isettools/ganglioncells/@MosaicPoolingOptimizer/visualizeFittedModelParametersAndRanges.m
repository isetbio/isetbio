function visualizeFittedModelParametersAndRanges(ax, modelParams, modelName)

    xTicks = 1:numel(modelParams.names);
    xTickLabels = modelParams.names;
    
    hold(ax, 'off');
    for iParam = 1:numel(modelParams.finalValues)
        
        plot(ax, iParam * [1 1], [0 1], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 4.0);
        hold(ax, 'on');

        if (modelParams.lowerBounds(iParam) == modelParams.upperBounds(iParam))
        else

            if (strcmp(modelParams.scaling, 'log'))
                fV = log10(modelParams.finalValues(iParam));
                fI = log10(modelParams.initialValues(iParam));
                lB = log10(modelParams.lowerBounds(iParam));
                uB = log10(modelParams.upperBounds(iParam));
                fValFinal = (fV - lB)/ (uB-lB);
                fValInitial = (fI - lB)/ (uB-lB);
            else
                fValFinal = (modelParams.finalValues(iParam) - modelParams.lowerBounds(iParam))/ ...
                            (modelParams.upperBounds(iParam)-modelParams.lowerBounds(iParam));
                fValInitial = (modelParams.initialValues(iParam) - modelParams.lowerBounds(iParam))/ ...
                            (modelParams.upperBounds(iParam)-modelParams.lowerBounds(iParam));
            end

            plot(ax, iParam * [1 1], [fValInitial fValFinal], 'b-', 'LineWidth', 3);
            plot(ax, iParam * [1 1], fValFinal, 'bo', 'MarkerSize', 14, 'MarkerFaceColor', [0.1 0.85 0.85], 'LineWidth', 1.5);
            text(ax, iParam+0.1, fValFinal, sprintf(' %1.3f', modelParams.finalValues(iParam)), 'FontSize', 13);
        end
        
        text(ax, iParam-0.25, -0.07, sprintf('%1.3f', modelParams.lowerBounds(iParam)), 'FontWeight', 'bold');
        text(ax, iParam-0.25, 1.05, sprintf('%1.3f', modelParams.upperBounds(iParam)), 'FontWeight', 'bold');
    end
    
    xtickangle(ax, 90);
    set(ax, 'XLim', [0.5 numel(modelParams.names)+0.7],  ...
            'XTick', xTicks, 'XTickLabel', xTickLabels, 'YLim', [-0.1 1.1], ...
            'YTickLabel', {}, 'FontSize', 14, 'YColor', 'none');
    
    axis (ax, 'square');
    ylabel(ax,'parameter value range');
    title(ax, modelName, 'FontSize', 16, 'Color', [0.3 0.3 0.3], 'FontWeight', 'normal');
    
end