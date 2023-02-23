function renderPerformance(ax, targetRsRcRatio, targetSCintSensRatio, ...
                        achievedRsRcRatio, achievedSCintSensRatio, ...
                        ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    targetKsKcRatio = targetSCintSensRatio / (targetRsRcRatio^2);
    achievedKsKcRatio = achievedSCintSensRatio/ (achievedRsRcRatio^2);

    bar(ax, 1, 100*((achievedRsRcRatio/targetRsRcRatio)-1), 'BaseValue', 0.0);
    hold(ax, 'on');
    bar(ax, 2, 100*((achievedSCintSensRatio/targetSCintSensRatio)-1), 'BaseValue', 0.0);
    bar(ax, 3, 100*((achievedKsKcRatio/targetKsKcRatio)-1), 'BaseValue', 0.0);
    
    grid(ax, 'on');
    set(ax, 'XLim', [0 4],  'XTick', [1 2 3], 'XTickLabel', {'Rs/Rc','S/C sens', 'Ks/Kc'}, ...
        'YLim', [-100 100], 'YTick', -100:20:100);
    ylabel(ax,'', 'FontAngle', ff.axisFontAngle);
    xlabel(ax, '', 'FontAngle', ff.axisFontAngle);

    xtickangle(ax, 90);
    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    % Plot title
    title(ax, sprintf('(achieved - target)/target'), ...
        'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
  
end
