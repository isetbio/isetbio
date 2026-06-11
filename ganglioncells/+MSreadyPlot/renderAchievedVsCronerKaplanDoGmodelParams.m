function renderAchievedVsCronerKaplanDoGmodelParams(ax, ...
    CK95dataX, CK95dataY, ISETBioDataX, ISETBioDataY, ...
    visualizedDataSetScatter, ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, YLims, YTicks, ...
    xAxisScaling, xAxisLabel, ...
    yAxisScaling, yAxisLabel, ...
    plotTitle, ff)

    cla(ax, 'reset');
    hold(ax, 'on');
    cellsNum = 0;

    for setIndex = 1:numel(ISETBioDataX)
        if (ismember(setIndex, visualizedDataSetScatter))
            cellsNum  = cellsNum  + numel(ISETBioDataX{setIndex}.data);
            scatter(ax, ISETBioDataX{setIndex}.data, ISETBioDataY{setIndex}.data, (ff.markerSize-8)^2,'o', ...
                'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', faceAlpha, 'MarkerEdgeAlpha', faceAlpha,  'LineWidth', ff.lineWidth);
        end

    end

    if (~isempty(CK95dataX)) && (~isempty(CK95dataY))
        % Add the Croner&Kaplan data
        scatter(ax, CK95dataX, ...
             CK95dataY, 144, 's', ...
             'MarkerFaceAlpha', 0.5, ...
             'MarkerEdgeColor', [0.2 0.2 0.2], ...
             'MarkerFaceColor', [0.8 0.8 0.8], ...
             'LineWidth', ff.lineWidth*0.75);
    end

    % Add the mean data
    
    
    if (strcmp(xAxisScaling, 'linear'))
        dX = (XLims(2)-XLims(1))*4/100;
    else
        dX = (XLims(2)-XLims(1))*0.1/100;
    end

    set(ax, 'XLim', [XLims(1) XLims(2)], 'YLim', [YLims(1) YLims(2)]);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'TickDir', 'both');
    set(ax, 'XScale', xAxisScaling, 'YScale', yAxisScaling);
    %axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);

    xlabel(ax, xAxisLabel, 'FontAngle', ff.axisFontAngle);
    ylabel(ax, yAxisLabel, 'FontAngle', ff.axisFontAngle);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    if (isempty(plotTitle))
        plotTitle = '';
    else
        plotTitle = sprintf('%d %s RGCs', cellsNum, plotTitle);
    end

    title(ax, plotTitle, 'fontSize', ff.titleFontSize, ...
         'Color', ff.titleColor, 'FontWeight', ff.titleFontWeight); 

end
