function [gridCoords,theRTVFobjIndicesForThisGrid] = subGridSpatialCoordsForConesNumPooled(theConesNumPooled, ...
    theConesNumPooledByTheRFcenterGrid, theOpticsPositionGrid)

    assert(ismember(theConesNumPooled, theConesNumPooledByTheRFcenterGrid), ...
        '%d center cone models are not computed for this mosaic', theConesNumPooled);
    
    theRTVFobjIndicesForThisGrid = find(theConesNumPooledByTheRFcenterGrid == theConesNumPooled);
    gridCoords = theOpticsPositionGrid(theRTVFobjIndicesForThisGrid,:);

end
