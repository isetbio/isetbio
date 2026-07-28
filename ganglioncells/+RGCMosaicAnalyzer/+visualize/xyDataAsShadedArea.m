function xyDataAsShadedArea(ax, x, y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
%
% RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
%
    x = [x fliplr(x)];
    y = [y y*0+baseline];

    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end