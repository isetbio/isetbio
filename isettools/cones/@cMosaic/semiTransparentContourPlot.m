% Function to generate a semitransparent controur plot
function semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor)
    % Compute contours at desired Z-level
    C = contourc(xSupport, ySupport, zData, zLevels);
    % Go through the contour matrix and plot each contour separately
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theNormalizedLevel = (theLevel-min(zLevels))/(max(zLevels)-min(zLevels));
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);

        lutIndex = 1+round(theNormalizedLevel*(size(cmap,1)-1));
        patch('Faces', f, 'Vertices', v, 'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(lutIndex,:), ...
            'FaceAlpha', alpha, ...
            'EdgeAlpha', alpha, ...
            'EdgeColor', contourLineColor, ...
            'LineStyle', '-', 'LineWidth', 1.0, ...
        'Parent', axesHandle);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end