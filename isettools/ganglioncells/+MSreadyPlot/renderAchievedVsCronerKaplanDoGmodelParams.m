function renderAchievedVsCronerKaplanDoGmodelParams(ax, ...
    CK95data, ISETBioData, visualizedISETBioDataSets, ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, YLims, YTicks, ...
    yAxisScaling, yAxisLabel, plotTitle, ff)

    cla(ax, 'reset');
    hold(ax, 'on');
    cellsNum = 0;

    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedISETBioDataSets))
            cellsNum  = cellsNum  + numel(ISETBioData{setIndex}.eccentricityDegs);
            scatter(ax, ISETBioData{setIndex}.eccentricityDegs, ISETBioData{setIndex}.data, (ff.markerSize-6)^2,'o', ...
                'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', ISETBioDataColor(setIndex,:), ...
                'MarkerFaceAlpha', faceAlpha, 'MarkerEdgeAlpha', 0.0,  'LineWidth', ff.lineWidth);
%             plot(ax, ISETBioData{setIndex}.eccentricityDegs, ISETBioData{setIndex}.data, '.', ...
%                  'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', ISETBioDataColor(setIndex,:), ...
%                  'MarkerSize', 6, 'LineWidth', 0.2);
        end

    end
    plot(ax, CK95data.eccentricityDegs, ...
             CK95data.data, 'ks', ...
             'MarkerSize', ff.markerSize, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', ff.lineWidth*0.75);
    
    
    dY = (YLims(2)-YLims(1))*4/100;
    
    set(ax, 'XLim', [XLims(1)-0.05 XLims(2)], 'YLim', [YLims(1)-dY YLims(2)]);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'TickDir', 'both');
    set(ax, 'XScale', 'log', 'YScale', yAxisScaling);
    axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);

    xlabel(ax,'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    ylabel(ax, yAxisLabel, 'FontAngle', ff.axisFontAngle);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    if (isempty(plotTitle))
        plotTitle = sprintf('n = %2.0f RGCs', cellsNum);
    end

    title(ax, plotTitle, 'fontSize', ff.titleFontSize, ...
         'Color', ff.titleColor, 'FontWeight', ff.titleFontWeight); 

end
