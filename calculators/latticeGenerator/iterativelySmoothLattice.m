function [rfPositions, rfPositionsHistory, maxMovements, iteration, terminationReason] = iterativelySmoothLattice(rfPositions, tabulatedSpacing, tabulatedEcc, iterativeParams, lambda, domain, visualizationParams)

    % Initiate state
    iteration = 0;
    lastTriangularizationIteration = 0;
    desiredSpringLengths = [];
    spacingDeviations = [];
    maxMovements = [];
    keepLooping = true;
    if (~isinf(iterativeParams.iterationsIntervalForSavingPositions))
        rfPositionsHistory(1,:,:) = single(rfPositions);
    else
        rfPositionsHistory = [];
    end
    
    tStart = clock;
    timePrevious = tStart;
    userRequestTerminationAtIteration = [];
    
    while (keepLooping)
        iteration = iteration + 1;
        
        % determine if we need to re-triangulate based on when we the
        % iterations that hapenned since the last one
        [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
            determineWhetherReTriangularizationIsNeeded(iteration, ...
            lastTriangularizationIteration, maxMovements, iterativeParams);
        
        % determine if we need to re-triangulate because local density is too high
        if (~reTriangulationIsNeeded)
            spacingDeviations = localRFSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc);
            reallyCloseRFsNum = numel(find(spacingDeviations >iterativeParams.thresholdSpacingDeviation));
    
            if (reallyCloseRFsNum > 0)
                reTriangulationIsNeeded = true;
                triangularizationTriggerEvent = 'rfs closer than thresholdSpacingDeviation';
            else
                reTriangulationIsNeeded = false;
                triangularizationTriggerEvent = '';
            end
        end
        
        if (reTriangulationIsNeeded)
            % Save iteration of triangularization
            lastTriangularizationIteration = iteration;
            
            % Re-triangulate and compute new spring structure
            [springs, springIndices, triangleIndices] = triangulate(rfPositions, lambda, domain);
        
            % Visualize lattice progression
            if (~visualizationParams.visualizeNothing)
                visualizeLatticeAndQuality(rfPositions, spacingDeviations, maxMovements, ...
                    triangleIndices, reTriangulationIsNeeded, triangularizationTriggerEvent, iteration, ...
                    iterativeParams, visualizationParams);
            end
        end % retriangularizationIsNeeded
        
        % Update rfPositions
        [rfPositions, desiredSpringLengths, maxMovements(iteration)] = ...
            updatePositions(rfPositions, desiredSpringLengths, springs, springIndices, ...
            tabulatedSpacing, tabulatedEcc, lambda, reTriangulationIsNeeded, domain, iterativeParams);

        % Save history        
        if ((mod(iteration,iterativeParams.iterationsIntervalForSavingPositions)==0) && ...
            (~isinf(iterativeParams.iterationsIntervalForSavingPositions)))
            rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
        end
        
        % Check different criteria for terminating looping
        if (maxMovements(iteration) < iterativeParams.dTolerance)
            keepLooping = false;
            terminationReason = 'movement less than specified tolerance';
            continue;
        end
        
        % Check whether we have achived the desired lattice qulity
        [keepLooping, minQualityValue] = checkForAdequateLatticeQuality(rfPositions, reTriangulationIsNeeded, triangleIndices, iterativeParams.minQValue);
        if (keepLooping == false)
            terminationReason = 'achieved specified lattice quality';
        end
        
        % Check whether we exceeded the max no of iterations
        if (iteration > iterativeParams.maxIterations)
            keepLooping = false;
            terminationReason = 'exceeded max no of iterations';
        end
        
        if (~isempty(userRequestTerminationAtIteration)) && (iteration >= userRequestTerminationAtIteration)
            keepLooping = false;
            terminationReason = sprintf('Terminated at iteration %d as per user request', iteration);
            continue;
        end
        
        if (visualizationParams.visualizeNothing)
            fprintf('Iteration %d (%2.2f hours): maxMovement = %2.4f microns, qVal = %2.3f\n', iteration, etime(clock, tStart)/60/60, maxMovements(iteration), minQualityValue);
        end
        
        
        % See if we need to query the user about terminating
        timeLapsedMinutes = etime(clock, timePrevious)/60;
        if (timeLapsedMinutes > iterativeParams.queryUserIntervalMinutes)
            fprintf('Another %d minute period has passed. Terminate soon?', iterativeParams.queryUserIntervalMinutes);
            userTermination = GetWithDefault(' If so enter # of iteration to terminate on. Otherwise hit enter to continue', 'continue');
            if (~strcmp(userTermination, 'continue'))
                userRequestTerminationAtIteration = str2double(userTermination);
                if (isnan(userRequestTerminationAtIteration))
                    userRequestTerminationAtIteration = [];
                end
            else
                fprintf('OK, will ask again in %d minutes.', iterativeParams.queryUserIntervalMinutes);
            end
            timePrevious = clock;
        end
            
    end % while keepLooping
    
    % Visualize final lattice progression
    if (~visualizationParams.visualizeNothing)
        triangleIndices = delaunayn(rfPositions);
        visualizeLatticeAndQuality(rfPositions, spacingDeviations, maxMovements, ...
            triangleIndices, reTriangulationIsNeeded, 'final iteration', iteration, ...
            iterativeParams, visualizationParams);
    end
            
end

function [keepLooping, minQualityValue] = checkForAdequateLatticeQuality(rfPositions, reTriangulationIsNeeded, triangleIndices, minQValue)
    
    % Compute quality values
    if (~reTriangulationIsNeeded)
        triangleIndices = delaunayn(rfPositions);
    end
    minQualityValue = computeHexLatticeQuality(rfPositions, triangleIndices);
    
    keepLooping = true;
    if (minQualityValue >= minQValue)
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


function [springs, springIndices, triangleIndices] = triangulate(rfPositions, lambda, domain)
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

function inputVal = GetWithDefault(prompt,defaultVal)
    if (ischar(defaultVal))
        inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
    else
        inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
    end
    if (isempty(inputVal))
        inputVal = defaultVal;
    end
end
