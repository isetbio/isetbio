function renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, ...
	edgeColor, faceColor, lineStyle, lineWidth)
% Render the patch array
%
% Syntax:
%   renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, ...
%       edgeColor, faceColor, lineStyle, lineWidth)
%
% Description:
%    Render the patch array
%
% Inputs:
%    axesHandle   - The handle to the patch axes
%    pixelOutline - The outline of the pixels
%    xCoords      - The X-axis coordinates
%    yCoords      - The Y-axis coordinates
%    edgeColor    - The color of the edge of the shape (line)
%    faceColor    - The color of the face of the shape (fill)
%    lineStyle    - The style of the lines (dash, solid, etc...)
%    lineWidth    - The width of the lines
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

verticesNum = numel(pixelOutline.x);
x = zeros(verticesNum, numel(xCoords));
y = zeros(verticesNum, numel(xCoords));

for vertexIndex = 1:verticesNum
    x(vertexIndex, :) = pixelOutline.x(vertexIndex) + xCoords;
    y(vertexIndex, :) = pixelOutline.y(vertexIndex) + yCoords;
end


patch(x, y, [0 0 0], 'EdgeColor', edgeColor, 'FaceAlpha', 1.0, ...
    'FaceColor', faceColor, 'LineWidth', lineWidth, ...
    'LineStyle', lineStyle, ...
    'Parent', axesHandle);

end