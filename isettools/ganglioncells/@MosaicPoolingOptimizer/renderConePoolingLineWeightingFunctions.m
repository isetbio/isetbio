function renderConePoolingLineWeightingFunctions(ax, ...
            centerLineWeightingFunction, surroundLineWeightingFunction, ...
            sensitivityRange)

    cMapEntries = 1024;
    cMap = (brewermap(cMapEntries, '*RdBu')).^1.0;


    centerLineWeightingFunction.spatialSupportDegs = (centerLineWeightingFunction.spatialSupportDegs(:))';
    centerLineWeightingFunction.amplitude = (centerLineWeightingFunction.amplitude(:))';
    surroundLineWeightingFunction.spatialSupportDegs = (surroundLineWeightingFunction.spatialSupportDegs(:))';
    surroundLineWeightingFunction.amplitude = (surroundLineWeightingFunction.amplitude(:))';
    
    xo = mean(centerLineWeightingFunction.spatialSupportDegs);
    x = centerLineWeightingFunction.spatialSupportDegs - xo;
    x1 = x(1);
    x2 = x(end);

    centerLineWeightingFunction.spatialSupportArcMin = (centerLineWeightingFunction.spatialSupportDegs  - xo)*60;
    surroundLineWeightingFunction.spatialSupportArcMin = (surroundLineWeightingFunction.spatialSupportDegs - xo)*60;

    shadedAreaBetweenTwoLines(ax, centerLineWeightingFunction.spatialSupportArcMin, ...
        centerLineWeightingFunction.amplitude, centerLineWeightingFunction.amplitude*0, ...
        cMap(512+256,:), 'none', 0.3, 1.5, '-');
    hold(ax, 'on');
    shadedAreaBetweenTwoLines(ax, surroundLineWeightingFunction.spatialSupportArcMin, ...
        -surroundLineWeightingFunction.amplitude, surroundLineWeightingFunction.amplitude*0, ...
        cMap(512-256,:), 'none', 0.3, 1.5, '-');

    plot(ax,surroundLineWeightingFunction.spatialSupportArcMin, ...
        centerLineWeightingFunction.amplitude - surroundLineWeightingFunction.amplitude, ...
        'k-', 'LineWidth', 1.0);
            
    plot(ax, surroundLineWeightingFunction.spatialSupportArcMin, surroundLineWeightingFunction.spatialSupportArcMin*0, 'k--', 'LineWidth', 0.5);
            
    xtickangle(ax, 0);
    set(ax, 'XTick', (-1: 0.05: 1)*60, 'XLim', [x1 x2]*60, 'YLim', sensitivityRange,  'TickDir', 'both');
    xlabel(ax, 'space, (arc min)');
    ylabel(ax, 'gain');
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'FontSize', 14)
end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
