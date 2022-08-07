function [insideBoundaryPointIndices, onBoundaryPointIndices] = pointsInsideBoundaryDefinedBySelectedPoints(...
    allPointPositions, selectedPointIndices)

    % Find indices of points are inside the boundary defined by a select subset of points

    allXcoords = squeeze(allPointPositions(:,1));
    allYcoords = squeeze(allPointPositions(:,2));
    
    shrinkFactor = 1.0;
    idx = boundary(allXcoords(selectedPointIndices), allYcoords(selectedPointIndices), shrinkFactor);

    boundingPolygonXcoords = allXcoords(selectedPointIndices(idx));
    boundingPolygonYcoords = allYcoords(selectedPointIndices(idx));
    [in,on] = inpolygon(allXcoords, allYcoords, boundingPolygonXcoords, boundingPolygonYcoords);

    insideBoundaryPointIndices = find((in == 1) | (on == 1));
    onBoundaryPointIndices = find(on == 1);
end