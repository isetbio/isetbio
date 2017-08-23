function renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
x = []; y = [];
triangleConeIndices = delaunayn([xHex(:), yHex(:)]);
for triangleIndex = 1:size(triangleConeIndices,1)
    coneIndices = triangleConeIndices(triangleIndex, :);
    xCoords = xHex(coneIndices);
    yCoords = yHex(coneIndices);
    for k = 1:numel(coneIndices)
        x = cat(2, x, xCoords);
        y = cat(2, y, yCoords);
    end
end
patch(x, y, [0 0 0], 'EdgeColor', meshEdgeColor, 'EdgeAlpha', meshEdgeAlpha, 'FaceAlpha', meshFaceAlpha, 'FaceColor', meshFaceColor, 'LineWidth', 1.5, 'LineStyle', lineStyle, 'Parent', axesHandle);
end