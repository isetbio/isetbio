function triangleIndices = triangularization(rfPositions, domainFunction, radius, borderTolerance)

    % Perform new Delaunay triangulation to determine the updated
    % topology of the truss.
    triangleIndices = delaunayn(rfPositions);

    % Compute the centroids of all triangles
    centroidPositions = 1.0/3.0 * (...
          rfPositions(triangleIndices(:, 1), :) + ...
          rfPositions(triangleIndices(:, 2), :) + ...
          rfPositions(triangleIndices(:, 3), :));

    % Remove centroids outside the desired region by applying the
    % signed distance function
    d = domainFunction(centroidPositions, radius);        
    triangleIndices = triangleIndices(d < borderTolerance, :);
    
end

