function renderSurroundMixHistograms(ax, lConeCenterSurroundConeMixRatios, ...
    mConeCenterSurroundConeMixRatios, ...
    plotTitle, ff, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTickLabel', false, @islogical);
    p.addParameter('noYTickLabel', false, @islogical);
    p.addParameter('targetRangeForSurroundConeMix', [], @isnumeric);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTickLabel = p.Results.noXTickLabel;
    noYTickLabel = p.Results.noYTickLabel;
    targetRangeForSurroundConeMix = p.Results.targetRangeForSurroundConeMix;
    
    % Compute histogram of surround cone mix ratios
    binWidth = 0.01;
    edges = 0:binWidth:(1.0+binWidth);
    LcenterCounts = histcounts(lConeCenterSurroundConeMixRatios, edges);
    McenterCounts = histcounts(mConeCenterSurroundConeMixRatios, edges);
    LMcounts = ([LcenterCounts(:) McenterCounts(:)]);
    LMcounts = LMcounts / sum(LMcounts(:));

    % The L-center RGCs
    [xx, yy1] = generateHistogramOutline(edges(1:size(LMcounts,1)), LMcounts(:,1));
    faceAlpha = 0.5;
    edgeColor = [0.5 0 0];
    faceColor = [1 0.35 0.35];
    shadedAreaBetweenTwoLines(ax, xx, yy1*100, yy1*0, ...
        faceColor, edgeColor, faceAlpha, 1.0, '-');
    hold(ax, 'on');

    % The M-center RGCs
    [xx, yy2] = generateHistogramOutline(edges(1:size(LMcounts,1)), LMcounts(:,2));
    edgeColor = [0 0.5 0];
    faceColor = [0.35 1.0 0.35];
    shadedAreaBetweenTwoLines(ax, xx, yy2*100, yy2*0, ...
        faceColor, edgeColor, faceAlpha, 1.0, '-');

    % The range of examined surround cone mix
    if (~isempty(targetRangeForSurroundConeMix))
        plot(ax, (targetRangeForSurroundConeMix(1)-0.25*binWidth)*[1 1], [0 100], 'k--', 'LineWidth', 1.5);
        plot(ax, (targetRangeForSurroundConeMix(2)+0.25*binWidth)*[1 1], [0 100], 'k--', 'LineWidth', 1.5);
    end

    grid(ax, 'on'); box(ax, 'off');
    Ymax = ceil(max([max(yy1) max(yy2)])*100);

    set(ax, 'XLim', [1*ff.axisOffsetFactor/2 1.0], 'YLim', [Ymax*ff.axisOffsetFactor/2 Ymax]);
    
    set(ax, 'XTick', 0:0.125:1, 'XTickLabel', {'0', '', '.25', '', '.5', '', '.75', '', '1'}, ...
            'YTick', 0:1:100, 'YTickLabel', sprintf('%2.0f\n', 0:0.5:100));
    
    set(ax, 'TickDir', 'both')

    % xy axis labels
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    if (~noXLabel)
        xlabel(ax, 'surround purity, \rho', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'percentage of cells (%)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', 0:0.1:1, 'YTickLabel', sprintf('%1.1f\n', 0:0.1:1));
    end

    if (~noYTickLabel)
    else
        set(ax, 'YTickLabel', {});
    end

    if (~noXTickLabel)
    else
        set(ax, 'XTickLabel', {});
    end


    % legends
    legend(ax, {' L-center: \rho = w^S_M/w^S_{LM} ', ' M-center:  \rho = w^S_L/w^S_{LM} '}, ...
        'Location', 'NorthOutside', 'NumColumns', 2, ...
        'FontSize', ff.legendFontSize, 'Box', 'off');

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

    drawnow;
end

function [xx, yy] = generateHistogramOutline(x,y)
    xx(1) = x(1);
    yy(1) = 0;

    for i = 1:(numel(x)-1)
        xx(numel(xx)+1) = x(i);
        yy(numel(yy)+1) = y(i);
        xx(numel(xx)+1) = x(i+1);
        yy(numel(yy)+1) = y(i);
    end

    xx(numel(xx)+1) = x(end);
    yy(numel(yy)+1) = 0;
end

 function p = shadedAreaBetweenTwoLines(ax,x, y1, y2, ...
     faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    y1 = reshape(y1, size(x));
    y2 = reshape(y2, size(x));

    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
 end