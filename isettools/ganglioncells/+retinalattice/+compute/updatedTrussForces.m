function [netForceVectors, desiredSpringLengths] = updatedTrussForces(rfPositions, tabulatedEcc, tabulatedRFspacing, springs, springIndices, desiredSpringLengths, rfSpacingFunctionFast, reTriangulationIsNeeded)
    % Compute spring vectors
    springVectors =  rfPositions(springs(:, 1), :) - rfPositions(springs(:, 2), :);
    % their centers
    springCenters = (rfPositions(springs(:, 1), :) + rfPositions(springs(:, 2), :)) / 2.0;
    % and their lengths
    springLengths = sqrt(sum(springVectors.^2, 2));
    
    if (reTriangulationIsNeeded)
        % Compute desired spring lengths. This is done by evaluating the
        % passed distance function at the spring centers.
        desiredSpringLengths = rfSpacingFunctionFast(springCenters, tabulatedEcc, tabulatedRFspacing);
    end
    
    % Normalize spring lengths
    normalizingFactor = sqrt(sum(springLengths .^ 2) / sum(desiredSpringLengths .^ 2));
    desiredSpringLengths = desiredSpringLengths * normalizingFactor;
   
    gain = 1.1;
    springForces = max(gain * desiredSpringLengths - springLengths, 0);

    % compute x, y-components of forces on each of the springs
    springForceXYcomponents = abs(springForces ./ springLengths * [1, 1] .* springVectors);

    % Compute net forces on each cone
    rfsNum = size(rfPositions,1);
    netForceVectors = zeros(rfsNum, 2);

    parfor rfIndex = 1:rfsNum
       % compute net force from all connected springs
       deltaPos = -bsxfun(@minus, springCenters(springIndices{rfIndex}, :), rfPositions(rfIndex, :));
       netForceVectors(rfIndex, :) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{rfIndex}, :), 1);
    end
        
end

