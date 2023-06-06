function renderAchievedVsCronerKaplanDoGmodelParamsHistograms(ax, ...
    CK95data, ISETBioData, visualizedISETBioDataSets, histogramEdges, ...
    ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, ...
    xAxisLabel, plotTitle, ff)

    
    cla(ax, 'reset');

    [counts,edges] = histcounts(CK95data.data, histogramEdges);
    countsPercentage = counts / numel(CK95data.data)*100;
    bar(ax, edges(1:end-1), countsPercentage, 1, ...
        'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none', ...
        'FaceAlpha', 1.0, 'LineWidth', ff.lineWidth);
    maxY = max(countsPercentage);

    hold(ax, 'on');
    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedISETBioDataSets))
            [counts,edges] = histcounts(ISETBioData{setIndex}.data, histogramEdges);
            countsPercentage = counts / numel(ISETBioData{setIndex}.data)*100;
            maxY = max([maxY max(countsPercentage)]);
            xOutline = [];
            yOutline = [];
            dd = (edges(2)-edges(1))/2;
            for i = 1:numel(edges)-1
                xOutline(numel(xOutline)+1) = edges(i)+dd;
                yOutline(numel(yOutline)+1) = countsPercentage(i);
                if (i < numel(edges)-1)
                xOutline(numel(xOutline)+1) = edges(i)+dd;
                yOutline(numel(yOutline)+1) = countsPercentage(i+1);
                end
            end

            bar(ax, edges(1:end-1), countsPercentage, 1, ...
                'FaceColor', ISETBioDataColor(setIndex,:), 'EdgeColor','none', ...
                'FaceAlpha', faceAlpha, 'LineWidth', ff.lineWidth);
            plot(ax, xOutline,yOutline, 'k-', 'LineWidth', 1.5, 'Color', ISETBioDataColor(setIndex,:)*0.5)
        end

    end

    yTickIncrement = 5;
    YLims(1) = 0;
    YLims(2) = ceil(maxY/yTickIncrement)*yTickIncrement;
    YTicks = 0:5:100;

    dX = (XLims(2)-XLims(1))*4/100;
    dY = (YLims(2)-YLims(1))*4/100;

    set(ax, 'XLim', [XLims(1)-dX XLims(2)], 'YLim', [YLims(1)-dY YLims(2)]);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'TickDir', 'both');
    set(ax, 'XScale', 'linear', 'YScale', 'linear');
    axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);

    xlabel(ax, xAxisLabel, 'FontAngle', ff.axisFontAngle);
    ylabel(ax, 'percentage', 'FontAngle', ff.axisFontAngle);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor,  'LineWidth', ff.axisLineWidth);
    
    title(ax, plotTitle, 'fontSize', ff.titleFontSize, ...
     'Color', ff.titleColor, 'FontWeight', ff.titleFontWeight); 

end
