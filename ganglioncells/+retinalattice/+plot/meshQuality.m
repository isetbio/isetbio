function subplotIndex = meshQuality(figNo,subplotIndex, histogramData, bin1Percent, iterationsHistory)
    if (isempty(figNo))
        hFig = figure(10); 
        subplotIndex = subplotIndex+1;
        if (subplotIndex >15)
            subplotIndex = 1;
        end
        if (subplotIndex == 1)
            clf;
            set(hFig, 'Position', [10 10 820 930]);
        end
        
        rows = 5; cols = 3;
        posVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rows, ...
           'colsNum', cols, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02);
        row = 1+floor((subplotIndex-1)/cols);
        col = 1+mod((subplotIndex-1),cols);
        subplot('Position', posVectors(row,col).v);
    end
 
    qLims = [0.6 1.005]; 
    bar(histogramData.x(2:end),histogramData.y,1); hold on;
    plot(bin1Percent(1)*[1 1], [0 max(histogramData.y)], 'r-', 'LineWidth', 1.5);
    plot(bin1Percent(end)*[1 1], [0 max(histogramData.y)], 'c-', 'LineWidth', 1.5);
    plot(bin1Percent(2)*[1 1], [0 max(histogramData.y)], 'k-',  'LineWidth', 1.5);
    plot(bin1Percent(3)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    plot(bin1Percent(4)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', qLims, 'YLim', [0 max(histogramData.y)], ...
        'XTick', [0.6:0.05:1.0],  'XTickLabel', {'.6', '', '.7', '', '.8', '', '.9', '', '1.'}, ...
        'YTickLabel', {}, 'FontSize', 12);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 12);
    if (isempty(figNo))
        title(sprintf('iteration:%d', iterationsHistory(end)))
        drawnow;
    end
    
    if (isempty(figNo))
        figure(11); hold on;
    end
    
end