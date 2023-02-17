function transparentContourPlot(axesHandle, xSupportMinutes, ySupportMinutes, zData, zLevels, cmap, varargin)

    p = inputParser;
    p.addParameter('faceAlpha', 0.4, @isscalar);
    p.parse(varargin{:});
    faceAlpha = p.Results.faceAlpha;

    C = contourc(xSupportMinutes, ySupportMinutes, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        patch(axesHandle, 'Faces', f, 'Vertices', v, 'EdgeColor', 0.5*(1-theLevel)*[1 1 1], ...
            'FaceColor', cmap(round(theLevel*size(cmap,1)),:), ...
            'FaceAlpha', faceAlpha, ...
            'LineStyle', '-', 'LineWidth', 1.0);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end