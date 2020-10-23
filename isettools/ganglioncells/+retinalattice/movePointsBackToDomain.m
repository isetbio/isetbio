function [rfPositions, d] = movePointsBackToDomain(rfPositions, lambdaMin, domainFunction, radius)
    % Find rfs that have moved outside the domain
    d = domainFunction(rfPositions, radius);
    outsideBoundaryIndices = d > 0;
        
    deps = sqrt(eps) * lambdaMin;
    
    % And project them back to the domain
    if (~isempty(outsideBoundaryIndices))
                % Compute numerical gradient along x-positions
                deltaX = [rfPositions(outsideBoundaryIndices, 1) + deps, rfPositions(outsideBoundaryIndices, 2)];
                dXgradient = (domainFunction(deltaX, radius) - d(outsideBoundaryIndices)) / deps;
            
                % Compute numerical gradient along y-positions
                deltaY = [rfPositions(outsideBoundaryIndices, 1), rfPositions(outsideBoundaryIndices, 2) + deps];
                dYgradient = (domainFunction(deltaY, radius) - d(outsideBoundaryIndices)) / deps;

                % Project these points back to boundary
                rfPositions(outsideBoundaryIndices, :) = rfPositions(outsideBoundaryIndices, :) - ...
                    [d(outsideBoundaryIndices) .* dXgradient, ...
                     d(outsideBoundaryIndices) .* dYgradient];
    end
        
end

