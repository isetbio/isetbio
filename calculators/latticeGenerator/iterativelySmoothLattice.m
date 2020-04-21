function [rfPositions] = iterativelySmoothLattice(rfPositions, tabulatedDensity, tabulatedEcc, iterativeParams, lambda, domain, neuronalType, whichEye)

    % Initiate state
    iteration = 0;
    lastTriangularizationIteration = 0;
    maxMovements = [];
    keepLooping = true;
    
    while (keepLooping)
        iteration = iteration + 1;
        
        % determine if we need to re-triangulate
        [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
            determineWhetherReTriangularizationIsNeeded(iteration, lastTriangularizationIteration, ...
            maxMovements, iterativeParams);
        
        if (reTriangulationIsNeeded)
            % Save iteration of triangularization
            lastTriangularizationIteration = iteration;
            
            % Re-triangulate and compute new spring structure
            [springs, springIndices] = triangulate(rfPositions, lambda, domain);
        end % retriangularizationIsNeeded
        
        % Update rfPositions
        [rfPositions, maxMovements(iteration)] = ...
            updatePositions(rfPositions, springs, springIndices, reTriangulationIsNeeded, domain, iterativeParams);
        
        % Check different criteria for terminating looping
        if (maxMovements(iteration) < iterativeParams.dTolerance)
            keepLooping = false; 
        end
        
        if (reTriangulationIsNeeded)
            [keepLooping, histogramData, minQualityValue] = ...
                checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions, triangleIndices, histogramWidths);
            visualizeHexLatticeQuality(histogramData, minQualityValue);
        end
        
         % Visualize lattice
        visualizeLattice(rfPositions);
        
    end % while keepLoopong
end

function [keepLooping, histogramData, minQualityValue] = ...
    checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions, triangleIndices)
    
    % Compute quality values
    [minQualityValue, histogramData] = computeHexLatticeQuality(rfPositions, triangleIndices);
    
    keepLooping = true;
    if (minQualityValue > 0.85)
        keepLooping = false;
    end
end


function [rfPositions, maxMovement] = updatePositions(rfPositions, springs, ...
    springIndices, reTriangulationIsNeeded, domain, iterativeParams)
        
    deltaT = 0.2;
    
    % Compute new spring vectors
    springVectors =  rfPositions(springs(:, 1), :) - rfPositions(springs(:, 2), :);
    % their centers
    springCenters = (rfPositions(springs(:, 1), :) + rfPositions(springs(:, 2), :)) / 2.0;
    % and their lengths
    springLengths = sqrt(sum(springVectors.^2, 2));

    if (reTriangulationIsNeeded)
        % Compute desired spring lengths. This is done by evaluating the
        % passed distance function at the spring centers.
        desiredSpringLengths = feval(gridParams.rfSpacingFunctionFast, springCenters, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons);
    end
        
    % Normalize spring lengths
    normalizingFactor = sqrt(sum(springLengths .^ 2) / ...
            sum(desiredSpringLengths .^ 2));
    desiredSpringLengths = desiredSpringLengths * normalizingFactor;
        
    gain = 1.1;
    springForces = max(gain * desiredSpringLengths - springLengths, 0);

    % compute x, y-components of forces on each of the springs
    springForceXYcomponents = abs(springForces ./ springLengths * [1, 1] .* springVectors);

    % Compute net forces on each cone
    netForceVectors = zeros(rfsNum, 2);
    parfor rfIndex = 1:rfsNum
       % compute net force from all connected springs
       deltaPos = -bsxfun(@minus, springCenters(springIndices{rfIndex}, :), rfPositions(rfIndex, :));
       netForceVectors(rfIndex, :) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{rfIndex}, :), 1);
    end
        
    % update cone positions according to netForceVectors
    rfPositions = rfPositions + deltaT * netForceVectors;
        
    % Project points that have gone outside the domain back inside
    rfPositions = projectPointsBackToEllipse(rfPositions, domain);
        
    % Check if all interior nodes move less than dTolerance
    movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d <-lambda/1000, :) .^2 , 2));
    maxMovement = prctile(movementAmplitudes, iterativeParams.maxMovementPercentile);
        
        
end

function rfPositions = projectPointsBackToEllipse(rfPositions, domain)
    % find RFs outside the domain
    d = feval(domain.function, rfPositions, domain.maxEcc, domain.ellipseAxes);
    idx = d > 0;
    
    % And project them back to the domain
    if (~isempty(outsideBoundaryIndices))
        % Compute numerical gradient along x-positions
        deltaX = [rfPositions(idx, 1)+deps, rfPositions(idx, 2)];
        deltaY = [rfPositions(idx, 1), rfPositions(idx, 2)+deps];
        dXgradient = 1/deps * (feval(domain.function, deltaX, domain.maxEcc, domain.ellipseAxes) - d(idx));
        dYgradient = 1/deps * (feval(domain.function, deltaY, domain.maxEcc, domain.ellipseAxes) - d(idx));

        % Project these points back to boundary
        rfPositions(idx, :) = rfPositions(idx, :) - [d(idx) .* dXgradient, d(idx) .* dYgradient];
    end
end


function [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
    determineWhetherReTriangularizationIsNeeded(iteration, lastTriangularizationIteration, ...
    maxMovements, iterativeParams)
 
    % We need to triangulate again if the movement in the current iteration was > the average movement in the last 2 iterations 
    if (numel(maxMovements)>3) && (maxMovements(iteration-1) > 0.52*(maxMovements(iteration-2)+maxMovements(iteration-3)))
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = 'movement stopped decreasing';
    end

    % We need to triangulate again if we went for maxIterationsToRetriangulate + some more since last triangularization
    if ((abs(lastTriangularizationIteration-iteration-1)) > iterativeParams.maxIterationsBeforeRetriangulation+min([10 round(iteration/20)]))
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = 'maxIterations passed';
    end

    % Do not triangulare if we did one less than minIterationsBeforeRetriangulation before
    if ((abs(lastTriangularizationAtIteration-iteration)) < iterativeParams.minIterationsBeforeRetriangulation+min([5 round(iteration/50)]))
        reTriangulationIsNeeded = false;
        triangularizationTriggerEvent = '';
    end

    %
    if (iteration==1)
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = '1st iteration';
    end
        
end


function [springs, springIndices] = triangulate(rfPositions, lambda, domain)
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
    d = feval(domain.function, centroidPositions, domain.maxEcc, domain.ellipseAxes);
    triangleIndices = triangleIndices(d < lambda/1000, :);

    % Create a list of the unique springs (each spring connecting 2 cones)
    springs = [...
            triangleIndices(:, [1, 2]); ...
            triangleIndices(:, [1, 3]); ...
            triangleIndices(:, [2, 3]) ...
    ];
    springs = unique(sort(springs, 2), 'rows');
           
    % find all springs connected to each rf
    rfsNum = size(rfPositions,1);
    springIndices = cell(1,rfsNum);
    for rfIndex = 1:rfsNum
        springIndices{rfIndex} = find((springs(:, 1) == rfIndex) | (springs(:, 2) == rfIndex));
    end
           
end
