function renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, ...
    meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
% Render the hex mesh for the cone mosaic hex object
%
% Syntax:
%    renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, ...
%        meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle)
%
% Description:
%    Render (draw) the hex mesh for the cone mosaic hex object using the
%    provided variables.
%
% Inputs:
%    axesHandle    - Handle to the cone mosaic hex axes
%    xHex          - x-bounded Hex
%    yHex          - y-bounded Hex
%    meshEdgeColor - The color of the edge of the mesh
%    meshFaceColor - The color for the face of the mesh
%    meshFaceAlpha - The alpha (transparency) value of the mesh faces
%    meshEdgeAlpha - The alpha value of the mesh edges
%    lineStyle     - The line style for the edging, ex. dash, solid, etc...
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
x = [];
y = [];
triangleConeIndices = delaunayn([xHex(:), yHex(:)]);
for triangleIndex = 1:size(triangleConeIndices, 1)
    coneIndices = triangleConeIndices(triangleIndex, :);
    xCoords = xHex(coneIndices);
    yCoords = yHex(coneIndices);
    for k = 1:numel(coneIndices)
        x = cat(2, x, xCoords);
        y = cat(2, y, yCoords);
    end
end
patch(x, y, [0 0 0], 'EdgeColor', meshEdgeColor, ...
    'EdgeAlpha', meshEdgeAlpha, 'FaceAlpha', meshFaceAlpha, ...
    'FaceColor', meshFaceColor, 'LineWidth', 1.0, ...
    'LineStyle', lineStyle, 'Parent', axesHandle);
end