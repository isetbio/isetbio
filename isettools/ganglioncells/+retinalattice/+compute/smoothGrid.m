function dataOut = smoothGrid(rfPositions, ...
        tabulatedEcc, tabulatedRFspacing, params, tStart)

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    
    % Params
    deltaT = 0.2;
    
    % State variables
    maxMovements = [];
    rfPositionsHistory = [];
    histogramWidths = [];
    reTriangulationIterations = [];
    
    % Iterative controls
    iteration = 0;
    lastTriangularizationIteration = 0;
    desiredSpringLengths = [];
    subplotIndex = 0;

    exceededMaxIterations = false;
    converged = false;
    terminateNowDueToReductionInLatticeQuality = false;
    userTerminated = false;
    
    while ((~exceededMaxIterations) && (~converged) && (~userTerminated) && ...
            (~terminateNowDueToReductionInLatticeQuality))
        
        iteration = iteration + 1;
        
        % Decide it triangulatization is needed
        [reTriangulationIsNeeded, triangularizationTriggerEvent] = ...
            retinalattice.check.triangulationConditions(maxMovements, ...
            iteration, lastTriangularizationIteration, ...
            params.minIterationsBeforeRetriangulation, ...
            params.maxIterationsBeforeRetriangulation);

        if (reTriangulationIsNeeded)
            % Save triangularization iteration
            lastTriangularizationIteration = iteration;
            
            % Save old rfPositions 
            oldRFPositions = rfPositions;
            
            % Compute new spring connections
            [springs, springIndices, triangleIndices] = retinalattice.compute.updatedTrussConnectivity(...
                rfPositions, params.domainFunction, params.radius, params.borderTolerance);
        end  % reTriangulationIsNeeded
        
        % Compute updated truss forces
        [netForceVectors, desiredSpringLengths] = retinalattice.compute.updatedTrussForces(...
            rfPositions, tabulatedEcc, tabulatedRFspacing, ...
            springs, springIndices, desiredSpringLengths, ...
            params.rfSpacingFastFunction, ...
            reTriangulationIsNeeded);
        
        % Compute updated rf positions using the updated truss forces
        rfPositions = rfPositions + deltaT * netForceVectors;
        
        % Move rfPositions back to the domain
        [rfPositions, d] = movePointsBackToDomain(...
            rfPositions, params.lambdaMinMicrons, params.domainFunction, params.radius);
        
        % Check if all nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d < -params.borderTolerance, :) .^2 , 2));
        maxMovements(iteration) = prctile(movementAmplitudes, params.maxMovementPercentile);

        % Check whether we converged due to dTolerance
        if (maxMovements(iteration) < params.dTolerance)
            converged = true; 
            terminationReason = sprintf('Max movement < dTolerance.');
        end
        
        % Check for early termination due to decrease in hex lattice quality
        if (reTriangulationIsNeeded)
        
            reTriangulationIterations = cat(2,reTriangulationIterations, iteration);
            
            % Compute mesh quality and see if quality is above desired
            % threshold and is not getting better
            
            [terminateNowDueToReductionInLatticeQuality, ...
             histogramData, histogramWidths, histogramDiffWidths, checkedBins] = ...
                retinalattice.check.earlyTerminationDueToHexLatticeQualityDecrease(...
                rfPositions, triangleIndices, histogramWidths, params.minHexQualityForTermination, params.qDistPercentile);
        
            if (isempty(rfPositionsHistory))
                rfPositionsHistory(1,:,:) = single(oldRFPositions);
                iterationsHistory = iteration;
            else
                rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(oldRFPositions), [1 size(oldRFPositions,1) size(oldRFPositions,2)]));
                iterationsHistory = cat(2, iterationsHistory, iteration);
            end
            
            retinalattice.plot.movementSequence([],maxMovements, params.dTolerance)
            subplotIndex = retinalattice.plot.meshQuality([],subplotIndex, histogramData, checkedBins, iterationsHistory);
        
            if (terminateNowDueToReductionInLatticeQuality)
                terminationReason = sprintf('Achieved min hex mesh quality.');
            end
            
        end
        
        if (mod(iteration,params.iterationIntervalForSavingPositions)==0)
            if (isempty(rfPositionsHistory))
                rfPositionsHistory(1,:,:) = single(rfPositions);
                iterationsHistory = iteration;
            else
                rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
                iterationsHistory = cat(2, iterationsHistory, iteration);
            end
        end
           
        if (iteration > params.maxIterations)
            exceededMaxIterations = true;
            terminationReason = sprintf('Exceeded max iterations (%d).', params.maxIterations);
        end
    end % while
    
    % Assemble data struct
    dataOut.rfPositions = rfPositions;
    dataOut.rfPositionsHistory = rfPositionsHistory;
    dataOut.iterationsHistory = iterationsHistory;
    dataOut.maxMovements = maxMovements;
    dataOut.reTriangulationIterations = reTriangulationIterations;
    dataOut.terminationReason = terminationReason;
end


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

