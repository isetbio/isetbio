function renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, edgeColor, faceColor, lineStyle)

verticesNum = numel(pixelOutline.x);
x = zeros(verticesNum, numel(xCoords));
y = zeros(verticesNum, numel(xCoords));

for vertexIndex = 1:verticesNum
    x(vertexIndex, :) = pixelOutline.x(vertexIndex) + xCoords;
    y(vertexIndex, :) = pixelOutline.y(vertexIndex) + yCoords;
end
patch(x, y, [0 0 0], 'EdgeColor', edgeColor, 'FaceAlpha', 1.0, 'FaceColor', faceColor, 'LineWidth', 1.0, 'LineStyle', lineStyle, 'Parent', axesHandle);
end