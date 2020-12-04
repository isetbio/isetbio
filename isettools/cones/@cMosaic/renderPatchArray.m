function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth)

    conesNum = numel(apertureRadii);
    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    if (numel(faceColors) == 1)
        colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
    else
        colors = [];
    end
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + rfCoords(coneIndex,2);
        if (numel(faceColors) == conesNum)
            colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
            size(colors)
        end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = 0.7;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
