function spacings = localRFSpacings(rfPositions, neighborsNum)
    % Find distances to neighors
    p = pdist2(rfPositions, rfPositions, 'euclidean', 'Smallest', neighborsNum+1);
    % exclude the points themselves which have 0 distance
    p = p(2:end,:);
    
    % Compute the median spacing among the neibhors
    spacings = median(p,1);
end
