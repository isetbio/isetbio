function renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(ax, ...
    CK95data, ISETBioData, visualizedDataSetScatter, visualizedDataSetMeans, ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, YLims, YTicks, ...
    yAxisScaling, yAxisLabel, plotTitle, ...
    employTemporalEquivalentEccentricity, ff)

    cla(ax, 'reset');
    hold(ax, 'on');
    cellsNum = 0;

    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedDataSetScatter))
            cellsNum  = cellsNum  + numel(ISETBioData{setIndex}.eccentricityDegs);
            p1 = scatter(ax, ISETBioData{setIndex}.eccentricityDegs, ISETBioData{setIndex}.data, (ff.markerSize-8)^2,'o', ...
                'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', faceAlpha, 'MarkerEdgeAlpha', faceAlpha,  'LineWidth', ff.lineWidth);
        end

    end

    % Add the Croner&Kaplan data
    p2 = scatter(ax, CK95data.eccentricityDegs, ...
             CK95data.data, 144, 's', ...
             'MarkerFaceAlpha', 0.5, ...
             'MarkerEdgeColor', [0.2 0.2 0.2], ...
             'MarkerFaceColor', [0.8 0.8 0.8], ...
             'LineWidth', ff.lineWidth*0.75);
    
    % Add the mean data
    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedDataSetMeans))
            X = ISETBioData{setIndex}.eccentricityDegs;
            Y = ISETBioData{setIndex}.data;
            % Compute mean values across eccentricities
            minEcc = 0.01;
            maxEcc = 30;
            [N,edges,bin] = histcounts(X, logspace(log10(minEcc), log10(maxEcc), 64));
            for iBin = 1:numel(edges)
                idx = find(bin == iBin);
                meanY(iBin) = mean(Y(idx));
            end

            saturatedColor = ISETBioDataColor(setIndex,:); %(ISETBioDataColor(setIndex,:)-min(ISETBioDataColor(setIndex,:)))./(max(ISETBioDataColor(setIndex,:))-min(ISETBioDataColor(setIndex,:)));

            plot(ax, edges, meanY, 'k-', 'LineWidth', 4, 'Color', [0 0 0]);
            plot(ax, edges, meanY, 'w--', 'LineWidth', 3, 'Color', saturatedColor);
        end
     end
    
    % Legend
    legend([p1 p2], {'ISETbio model', 'Croner&Kaplan ''95'}, 'Location', 'NorthWest', 'box', 'off', 'Color', 'none');

    dY = (YLims(2)-YLims(1))*4/100;
    
    set(ax, 'XLim', [XLims(1)-0.05 XLims(2)], 'YLim', [YLims(1)-dY YLims(2)]);
    set(ax, 'XTick', XTicks, 'YTick', YTicks);
    set(ax, 'TickDir', 'both');
    set(ax, 'XScale', 'log', 'YScale', yAxisScaling);
    %axis(ax, 'square');
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