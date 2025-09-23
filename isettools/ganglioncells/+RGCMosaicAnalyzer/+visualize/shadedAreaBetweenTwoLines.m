function p = shadedAreaBetweenTwoLines(ax,x, y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
%
% RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
%
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