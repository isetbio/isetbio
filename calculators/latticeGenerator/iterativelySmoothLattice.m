function [rfPositions] = iterativelySmoothLattice(rfPositions, tabulatedSpacing, tabulatedEcc, iterativeParams, lambda, domain)

    % Initiate state
    iteration = 0;
    lastTriangularizationIteration = 0;
    desiredSpringLengths = [];
    spacingDeviations = [];
    maxMovements = [];
    keepLooping = true;
    visualizeLatticeGridQuality = true;
    
    if (visualizeLatticeGridQuality)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', 28, ...
                'figureHeightInches', 14);
    end
    
    while (keepLooping)
        iteration = iteration + 1;
        
        % determine if we need to re-triangulate based on when we the
        % iterations that hapenned since the last one
        [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
            determineWhetherReTriangularizationIsNeeded(iteration, ...
            lastTriangularizationIteration, maxMovements, iterativeParams);
        
        % determine if we need to re-triangulate because local density is
        % too high
        if (~reTriangulationIsNeeded)
            [reTriangulationIsNeeded, triangularizationTriggerEvent, spacingDeviations] = ...
                checkForLocalSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc, iterativeParams);
        end
        
        if (reTriangulationIsNeeded)
            % Save iteration of triangularization
            lastTriangularizationIteration = iteration;
            
            % Re-triangulate and compute new spring structure
            [springs, springIndices] = triangulate(rfPositions, lambda, domain);
        end % retriangularizationIsNeeded
        
        if (visualizeLatticeGridQuality)
            visualizeLatticeAndQuality(rfPositions, spacingDeviations, reTriangulationIsNeeded, iteration);
        end
        
        % Update rfPositions
        [rfPositions, desiredSpringLengths, maxMovements(iteration)] = ...
            updatePositions(rfPositions, desiredSpringLengths, springs, springIndices, ...
            tabulatedSpacing, tabulatedEcc, lambda, reTriangulationIsNeeded, domain, iterativeParams);

        % Check different criteria for terminating looping
        if (maxMovements(iteration) < iterativeParams.dTolerance)
            keepLooping = false; 
        end
        
        if (reTriangulationIsNeeded)
            [keepLooping, histogramData, minQualityValue] = ...
                checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions);
            visualizeHexLatticeQuality(histogramData, minQualityValue);
        end
        
        figure(55);
        plot(1:iteration, maxMovements,'ks-');
        drawnow;
        
        % Visualize lattice
        %visualizeLattice(rfPositions);
        
        if (iteration > iterativeParams.maxIterations)
            keepLooping = false;
        end
    end % while keepLoopong
end

function [keepLooping, histogramData, minQualityValue] = ...
    checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions)
    
    % Compute quality values
    triangleIndices = delaunayn(rfPositions);
    [minQualityValue, histogramData] = computeHexLatticeQuality(rfPositions, triangleIndices);
    
    keepLooping = true;
    if (minQualityValue > 0.8)
        keepLooping = false;
    end
end


function [rfPositions, desiredSpringLengths, maxMovement] = updatePositions(rfPositions, desiredSpringLengths, springs, ...
    springIndices, tabulatedSpacing, tabulatedEcc, lambda, reTriangulationIsNeeded, domain, iterativeParams)
        
    deltaT = 0.2;
    rfsNum = size(rfPositions,1);
    
    % Compute new spring vectors
    springVectors =  rfPositions(springs(:, 1), :) - rfPositions(springs(:, 2), :);
    % their centers
    springCenters = (rfPositions(springs(:, 1), :) + rfPositions(springs(:, 2), :)) / 2.0;
    % and their lengths
    springLengths = sqrt(sum(springVectors.^2, 2));

    if (reTriangulationIsNeeded)
        % Compute desired spring lengths. 
        neighborsNum = 9;
        desiredSpringLengths = lookUpValues(springCenters, tabulatedEcc, tabulatedSpacing, neighborsNum);
    end
        

    % Normalize spring lengths
    normalizingFactor = sqrt(sum(springLengths .^ 2) / sum(desiredSpringLengths .^ 2));
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
    [rfPositions, d] = projectPointsBackToEllipse(rfPositions, lambda, domain);
        
    % Compute max movement of all interior nodes
    movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d <-lambda/1000, :) .^2 , 2));
    maxMovement = prctile(movementAmplitudes, iterativeParams.maxMovementPercentile);   
end

function [rfPositions, d] = projectPointsBackToEllipse(rfPositions, lambda, domain)
    % find RFs outside the domain
    d = feval(domain.function, rfPositions, domain.maxEcc, domain.ellipseAxes);
    idx = d > 0;
    deps = sqrt(eps) * lambda;
    
    % And project them back to the domain
    if (~isempty(idx))
        % Compute numerical gradient along x-positions
        deltaX = [rfPositions(idx, 1)+deps, rfPositions(idx, 2)];
        deltaY = [rfPositions(idx, 1), rfPositions(idx, 2)+deps];
        dXgradient = 1/deps * (feval(domain.function, deltaX, domain.maxEcc, domain.ellipseAxes) - d(idx));
        dYgradient = 1/deps * (feval(domain.function, deltaY, domain.maxEcc, domain.ellipseAxes) - d(idx));

        % Project these points back to boundary
        rfPositions(idx, :) = rfPositions(idx, :) - [d(idx) .* dXgradient, d(idx) .* dYgradient];
    end
end

function [reTriangulationIsNeeded, triangularizationTriggerEvent, spacingDeviations] = checkForLocalSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc, iterativeParams)

    % Find distances to neighors
    neighborsNum = 1;
    spacings = localRFSpacings(rfPositions, neighborsNum);
    desiredSpacings = (lookUpValues(rfPositions, tabulatedEcc, tabulatedSpacing, neighborsNum))';
    spacingDeviations = (abs(desiredSpacings-spacings))./desiredSpacings;
    reallyCloseRFsNum = numel(find(spacingDeviations >iterativeParams.thresholdSpacingDeviation));
    
    if (reallyCloseRFsNum > 0)
        
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = 'rfs closer than thresholdSpacingDeviation';
        
        visualizeDeviationMap = true;
        if (visualizeDeviationMap)
            % Sampling vector
            sampling = struct('minPos', 1, 'maxPos', max(abs(rfPositions(:))), 'intervals', 100, 'scale', 'log');
            % Generate 2D map from scattered values
            [deviationMap, mapSupport] = mapFromScatteredPositions(rfPositions, spacingDeviations, sampling);
            
            figure(444);
            contourf(mapSupport(:,:,1), mapSupport(:,:,2), deviationMap, 0:0.05:1.0);
            set(gca, 'CLim', [0 1], 'ZLim', [0 1]);
            colormap(jet)
            axis 'square'
            colorbar
            title(sprintf('deviations = %f-%f', min(deviationMap(:)), max(deviationMap(:))));
        end
        
    
    else
        reTriangulationIsNeeded = false;
        triangularizationTriggerEvent = '';
    end
end


function [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
    determineWhetherReTriangularizationIsNeeded(iteration, lastTriangularizationIteration, ...
    maxMovements, iterativeParams)
 
    % Start with no re-triangulatization
    reTriangulationIsNeeded = false;
    triangularizationTriggerEvent = '';
    
    if (iteration==1)
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = '1st iteration';
        return;
    end
    
    
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
    if ((abs(lastTriangularizationIteration-iteration)) < iterativeParams.minIterationsBeforeRetriangulation+min([5 round(iteration/50)]))
        reTriangulationIsNeeded = false;
        triangularizationTriggerEvent = '';
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
