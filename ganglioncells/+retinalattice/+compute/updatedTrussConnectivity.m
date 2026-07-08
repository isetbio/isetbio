function [springs, springIndices, triangleIndices] = updatedTrussConnectivity(...
    rfPositions, domainFunction, radius, borderTolerance)
    
    triangleIndices = retinalattice.compute.triangularization(rfPositions, ...
        domainFunction, radius, borderTolerance);
    
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

