function idx = indicesOfPointsOutside(obj, points)

    assert(size(points,2) == 2, 'points must be an [Nx2] matrix');

    idx = setdiff(1:size(points,1), obj.indicesOfPointsInside(points));
end
