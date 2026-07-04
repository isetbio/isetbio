function renderPerformance(ax, targetRsRcRatio, targetSCintSensRatio, ...
                        achievedRsRcRatio, achievedSCintSensRatio, ...
                        ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noTitle', false, @islogical);

    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noTitle = p.Results.noTitle;

    targetKsKcRatio = targetSCintSensRatio / (targetRsRcRatio^2);
    achievedKsKcRatio = achievedSCintSensRatio/ (achievedRsRcRatio^2);

    bar(ax, 1, 100*((achievedRsRcRatio/targetRsRcRatio)-1), 0.5, 'BaseValue', 0.0);
    hold(ax, 'on');
    bar(ax, 2, 100*((achievedSCintSensRatio/targetSCintSensRatio)-1), 0.5, 'BaseValue', 0.0);
    bar(ax, 3, 100*((achievedKsKcRatio/targetKsKcRatio)-1), 0.5, 'BaseValue', 0.0);
    
    grid(ax, 'on');
    box(ax, 'off');
    
    if (~noYLabel)
        ylabel(ax,'residual (%)', 'FontAngle', ff.axisFontAngle);
    end

    xlabel(ax, '', 'FontAngle', ff.axisFontAngle);

    
    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    % 0-deg roated x-ticks
    xtickangle(ax, 0);

    set(ax, 'XLim', [0.5+3*ff.axisOffsetFactor 3.5],  'XTick', 0.5:0.5:3.5, ...
        'XTickLabel', {'', 'R_s/R_c', '', 'intS_s/intS_c', '', 'K_s/K_c', ''}, ...
        'YLim', [-80+200*ff.axisOffsetFactor  80], ...
        'YTick', -100:20:100);

    % Plot title
    if (~noTitle)
        title(ax, sprintf('(model - Croner&Kaplan)/Croner&Kaplan'), ...
            'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end
