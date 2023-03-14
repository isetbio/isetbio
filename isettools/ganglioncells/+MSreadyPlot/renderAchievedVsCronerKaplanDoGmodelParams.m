function renderAchievedVsCronerKaplanDoGmodelParams(ax, ...
    CK95data, ISETBioData, visualizedDataSet, ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, YLims, YTicks, ...
    yAxisScaling, yAxisLabel, plotTitle, ...
    employTemporalEquivalentEccentricity, ff)

    cla(ax, 'reset');
    hold(ax, 'on');
    cellsNum = 0;

    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedDataSet))
            cellsNum  = cellsNum  + numel(ISETBioData{setIndex}.eccentricityDegs);
            scatter(ax, ISETBioData{setIndex}.eccentricityDegs, ISETBioData{setIndex}.data, (ff.markerSize-6)^2,'o', ...
                'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', ISETBioDataColor(setIndex,:), ...
                'MarkerFaceAlpha', faceAlpha, 'MarkerEdgeAlpha', 0.0,  'LineWidth', ff.lineWidth);
        end

    end

    % Add the Croner&Kaplan data
    scatter(ax, CK95data.eccentricityDegs, ...
             CK95data.data, 144, 's', ...
             'MarkerFaceAlpha', 0.5, ...
             'MarkerEdgeColor', [0.2 0.2 0.2], ...
             'MarkerFaceColor', [0.8 0.8 0.8], ...
             'LineWidth', ff.lineWidth*0.75);
    
    % Add the mean data
    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedDataSet))
            X = ISETBioData{setIndex}.eccentricityDegs;
            Y = ISETBioData{setIndex}.data;
            % Compute mean values across eccentricities
            minEcc = 0.3;
            maxEcc = 30;
            [N,edges,bin] = histcounts(X, logspace(log10(minEcc), log10(maxEcc), 32));
            for iBin = 1:numel(edges)
                idx = find(bin == iBin);
                meanY(iBin) = mean(Y(idx));
            end

            saturatedColor = (ISETBioDataColor(setIndex,:)-min(ISETBioDataColor(setIndex,:)))./(max(ISETBioDataColor(setIndex,:))-min(ISETBioDataColor(setIndex,:)));

            plot(ax, edges, meanY, 'k-', 'LineWidth', 4, 'Color', ISETBioDataColor(setIndex,:));
            plot(ax, edges, meanY, 'w-', 'LineWidth', 2, 'Color', saturatedColor);
        end
     end
    
    
    dY = (YLims(2)-YLims(1))*4/100;
    
    set(ax, 'XLim', [XLims(1)-0.05 XLims(2)], 'YLim', [YLims(1)-dY YLims(2)]);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'TickDir', 'both');
    set(ax, 'XScale', 'log', 'YScale', yAxisScaling);
    axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);

    if (employTemporalEquivalentEccentricity)
        xlabel(ax,'temporal equivalent eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    else
        xlabel(ax,'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    end

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
