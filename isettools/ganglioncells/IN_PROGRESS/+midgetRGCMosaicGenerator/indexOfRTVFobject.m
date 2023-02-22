function theRTVFobjectIndex = indexOfRTVFobject(...
    conesNum, positionDegs, ...
    theConesNumPooledByTheRFcenterGrid, ...
    theOpticsPositionGrid)

    centerConeMatchObjIndices = find(theConesNumPooledByTheRFcenterGrid == conesNum);
    distancesToSamplingGridPositions = sqrt(sum((bsxfun(@minus, theOpticsPositionGrid(centerConeMatchObjIndices,:), positionDegs)).^2,2));
    [~, idx] = min(distancesToSamplingGridPositions);
    theRTVFobjectIndex = centerConeMatchObjIndices(idx);

end