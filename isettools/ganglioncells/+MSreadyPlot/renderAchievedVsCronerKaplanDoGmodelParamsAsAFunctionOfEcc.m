function hFig = renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(ax, ...
    CK95data, ISETBioData, visualizedDataSetScatter, visualizedDataSetMeans, ISETBioDataColor, faceAlpha, ...
    XLims, XTicks, YLims, YTicks, ...
    yAxisScaling, yAxisLabel, plotTitle, ...
    showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
    employTemporalEquivalentEccentricity, ff)


    minEcc = 0.01;
    maxEcc = 30;
    eccBins = 64;
    eccSupportForISETBioMeanData = logspace(log10(minEcc), log10(maxEcc), eccBins);

    % Ecc supports
    minEccForZscore = max([min(eccSupportForISETBioMeanData) min(CK95data.eccentricityDegs)]);
    maxEccForZscore = min([max(eccSupportForISETBioMeanData) max(CK95data.eccentricityDegs)]);
    eccBinsForZscore = 10;
    eccSupportForZscores = logspace(log10(minEccForZscore), log10(maxEccForZscore), eccBinsForZscore);


    if (isempty(ax))
        hFig = figure(1); clf;
        if (showZscoresInsteadOfData)
            ff = MSreadyPlot.figureFormat('1x1 small');
        else
            ff = MSreadyPlot.figureFormat('1x1 small tall');
        end
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);
        ax = theAxes{1,1};
    end   

    cla(ax, 'reset');
    hold(ax, 'on');
    cellsNum = 0;

    for setIndex = 1:numel(ISETBioData)
        if (ismember(setIndex, visualizedDataSetScatter))
            cellsNum  = cellsNum  + numel(ISETBioData{setIndex}.eccentricityDegs);
            if (~showZscoresInsteadOfData) && (~onlyShowCronerKaplan95Data)
                p1 = scatter(ax, ISETBioData{setIndex}.eccentricityDegs, ISETBioData{setIndex}.data, (ff.markerSize-8)^2,'o', ...
                    'MarkerFaceColor', ISETBioDataColor(setIndex,:), 'MarkerEdgeColor', 'none', ...
                    'MarkerFaceAlpha', faceAlpha, 'MarkerEdgeAlpha', faceAlpha,  'LineWidth', ff.lineWidth);
            end
        end
        
        
        if (showZscoresInsteadOfData)
            % Compute z-scores
            zScores = cell(1,numel(eccSupportForZscores));
            ZscoreEccentricity = nan(1,numel(eccSupportForZscores));

            for iEcc = 1:numel(eccSupportForZscores)
                if (iEcc>1)
                    eccLeft = eccSupportForZscores(iEcc-1);
                else
                    eccLeft = eccSupportForZscores(1)-0.1;
                end
        
                if (iEcc<numel(eccSupportForZscores))
                    eccRight = eccSupportForZscores(iEcc+1);
                else
                    eccRight = eccSupportForZscores(end)*1.1;
                end
        
                idx = find((CK95data.eccentricityDegs>=eccLeft)&(CK95data.eccentricityDegs<=eccRight));
                if (numel(idx)>1)
                    theCK95values = CK95data.data(idx);
                    meanOfCK95values = mean(theCK95values);
                    standardDeviationOfCK95values = std(theCK95values);
        
                    idx = find((ISETBioData{setIndex}.eccentricityDegs>=eccLeft)&(ISETBioData{setIndex}.eccentricityDegs<=eccRight));
                    theISETBiovalues = ISETBioData{setIndex}.data(idx);
                    allZscores = (theISETBiovalues-meanOfCK95values)/standardDeviationOfCK95values;
                    zScores{iEcc} = struct('mean', mean(allZscores), 'prctiles', prctile(allZscores(:),[1 99]));
                    ZscoreEccentricity(iEcc) = mean([eccLeft eccRight]);
                end
            end
            maxZscore = 2.5;
        end
    end

    if (showZscoresInsteadOfData)
        for iEcc = 1:numel(eccSupportForZscores)
            hold(ax, 'on');
            if (~isempty(zScores{iEcc}))
                mm = zScores{iEcc}.mean;
                zScoreRange = zScores{iEcc}.prctiles;
                plot(ax, ZscoreEccentricity(iEcc)*[1 1], zScoreRange, 'r-', ...
                     'LineWidth', ff.lineWidth);
                plot(ax, ZscoreEccentricity(iEcc)*[1 1], mm * [1 1], 'ro', ...
                    'MarkerSize', ff.markerSize, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', ff.lineWidth);
                
            end
        end

        % Font size
        set(ax, 'FontSize', ff.fontSize);
    
        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
        

        set(ax, 'XLim', [XLims(1)-0.05 XLims(2)], 'YLim', maxZscore*[-1.0 1.0], 'YTick', -3:0.5:3, 'YTickLabel', sprintf('%2.1f\n', -3:0.5:3));
        set(ax, 'XTick', XTicks);
        set(ax, 'TickDir', 'both');
        set(ax, 'XScale', 'log', 'YScale', 'linear');

        grid(ax, 'on');
        box(ax, 'off');
        xtickangle(ax, 0);

        if (employTemporalEquivalentEccentricity)
            xlabel(ax,'temporal equivalent eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
        else
            xlabel(ax,'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
        end
        
        ylabel(ax,'z-score', 'FontAngle', ff.axisFontAngle);

    else

        % Add the Croner&Kaplan data
        p2 = scatter(ax, CK95data.eccentricityDegs, ...
                 CK95data.data, ff.markerSize^2, 's', ...
                 'MarkerFaceAlpha', 0.5, ...
                 'MarkerEdgeColor', [0.2 0.2 0.2], ...
                 'MarkerFaceColor', [0.8 0.8 0.8], ...
                 'LineWidth', ff.lineWidth*0.75);
        
        if (~onlyShowCronerKaplan95Data)
            % Add the mean ISETBio data
            for setIndex = 1:numel(ISETBioData)
                if (ismember(setIndex, visualizedDataSetMeans))
                    X = ISETBioData{setIndex}.eccentricityDegs;
                    Y = ISETBioData{setIndex}.data;
                    % Compute mean values across eccentricities
                    [N,edges,bin] = histcounts(X, eccSupportForISETBioMeanData);
                    for iBin = 1:numel(edges)
                        idx = find(bin == iBin);
                        meanY(iBin) = mean(Y(idx));
                    end
                    plot(ax, edges, meanY, 'k-', 'LineWidth', 4, 'Color', [0 0 0]);
                    plot(ax, edges, meanY, 'w--', 'LineWidth', 3, 'Color', ISETBioDataColor(setIndex,:));
                end
            end
        end

        
        % Legend
        if (onlyShowCronerKaplan95Data)
            legend(p2, {'Croner&Kaplan ''95'}, 'Location', 'NorthWest', 'box', 'off', 'Color', 'none');
        else
            legend([p1 p2], {'ISETbio model', 'Croner&Kaplan ''95'}, 'Location', 'NorthWest', 'box', 'off', 'Color', 'none');
        end

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

end
