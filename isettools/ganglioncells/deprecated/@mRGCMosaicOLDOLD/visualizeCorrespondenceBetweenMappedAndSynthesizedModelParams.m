function visualizeCorrespondenceBetweenMappedAndSynthesizedModelParams(obj,visuallyMappedModelParams)

    % Retrieve the visual synthesized params
    synthesizedParams = obj. synthesizedRFparams.visual;
    
    % Compute synthesized ratios
    synthesizedSurroundToCenterRadiusRatios = ...
        synthesizedParams.surroundCharacteristicRadiiDegs./synthesizedParams.centerCharacteristicRadiiDegs;
    
    synthesizedSurroundToCenterPeakSensitivityRatios = ...
        synthesizedParams.surroundPeakSensitivities./synthesizedParams.centerPeakSensitivities;
    
    synthesizedSurroundToCenterIntegratedSensRatios = ...
        synthesizedSurroundToCenterRadiusRatios.^2 .* synthesizedSurroundToCenterPeakSensitivityRatios;
    
    % Compute visually-mapped ratios
    visuallyMappedSurroundToCenterRadiusRatios = ...
        visuallyMappedModelParams.surroundCharacteristicRadiiDegs./visuallyMappedModelParams.centerCharacteristicRadiiDegs;
    
    visuallyMappedSurroundToCenterPeakSensitivityRatios = ...
        visuallyMappedModelParams.surroundPeakSensitivities./visuallyMappedModelParams.centerPeakSensitivities;
    
    visuallyMappedSurroundToCenterIntegratedSensRatios = ...
        visuallyMappedSurroundToCenterRadiusRatios.^2 .* visuallyMappedSurroundToCenterPeakSensitivityRatios;
    
    
    % Plot everything
    figure(1111); clf;
    ax = subplot(2,3,1);
    correspondencePlot(ax, synthesizedParams.centerCharacteristicRadiiDegs*60, ...
         visuallyMappedModelParams.centerCharacteristicRadiiDegs*60, [0 4],...
         'linear', 'Rc (arc min)');
    
    ax = subplot(2,3,2);
    correspondencePlot(ax, synthesizedParams.surroundCharacteristicRadiiDegs*60, ...
         visuallyMappedModelParams.surroundCharacteristicRadiiDegs*60, [0 30],...
         'linear', 'Rs (arc min)');
    

    ax = subplot(2,3,3);
    correspondencePlot(ax, synthesizedSurroundToCenterRadiusRatios, ...
         visuallyMappedSurroundToCenterRadiusRatios, [0 10],...
         'linear', 'Rs/Rc');

    
    ax = subplot(2,3,4);
    correspondencePlot(ax, synthesizedSurroundToCenterPeakSensitivityRatios, ...
         visuallyMappedSurroundToCenterPeakSensitivityRatios, [0.001 0.1],...
         'log', 'Ks/Kc');
    
    ax = subplot(2,3,5);
    correspondencePlot(ax, synthesizedSurroundToCenterIntegratedSensRatios, ....
        visuallyMappedSurroundToCenterIntegratedSensRatios, [0 1], ...
        'linear', 'integrated S/C sensitivity');
end

function correspondencePlot(ax, x, y, xyRange, axisScaling, variableName)
    plot(ax, [xyRange(1) xyRange(2)], [xyRange(1) xyRange(2)], 'k-', 'LineWidth', 1.5); hold(ax, 'on');
    plot(ax, [0 xyRange(2)], [0 xyRange(2)*2], 'k--');
    plot(ax, [0 xyRange(2)], [0 xyRange(2)*0.5], 'k--');
    [N,edges] = histcounts(y);
    barh(ax,edges(1:end-1),N/max(N)*(xyRange(2)-xyRange(1)));
    plot(ax, x, y,'ro', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
    hold(ax, 'off');
    axis(ax, 'square');
    set(ax, 'XLim', xyRange, 'YLim', xyRange, 'XScale', axisScaling, 'YScale', axisScaling);
    xlabel(ax,'synthesized');
    ylabel(ax,'mapped');
    title(ax,variableName);
    
end

