function [springs, springIndices, triangleIndices] = updateTrussConnectivity(...
    rfPositions, domainFunction, radius, borderTolerance)
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

    % Create a list of the unique springs (each spring connecting 2 cones)
    springs = [...
            triangleIndices(:, [1, 2]); ...
            triangleIndices(:, [1, 3]); ...
            triangleIndices(:, [2, 3]) ...
    ];
    springs = unique(sort(springs, 2), 'rows');
           
    % find all springs connected to this cone
    rfsNum = size(rfPositions,1);
    springIndices = cell(1,rfsNum);
    for rfIndex = 1:rfsNum
       springIndices{rfIndex} = find((springs(:, 1) == rfIndex) | (springs(:, 2) == rfIndex));
    end
           
end

