function [reTriangulationIsNeeded, triangularizationTriggerEvent] = triangulationConditions(...
    maxMovements, iteration, lastTriangularizationIteration, ...
    minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation)

    reTriangulationIsNeeded = false;
    triangularizationTriggerEvent = '';
    
    % Triangulate again if the movement in the current iteration was > the average movement in the last 2 iterations 
    if (numel(maxMovements)>3) && (maxMovements(iteration-1) > 1.02*0.5*(maxMovements(iteration-2)+maxMovements(iteration-3)))
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = 'lots of node movement';
    end
        
    % Triangulate again if we went for maxIterationsToRetriangulate + some more since last triangularization
    if ((abs(lastTriangularizationIteration-iteration-1)) > maxIterationsBeforeRetriangulation+min([10 round(iteration/20)]))
        reTriangulationIsNeeded = true;
        triangularizationTriggerEvent = 'exceeded max iterations between triangularizations';
    end
        
    % Do not triangulare again if we triangulated recently
    if ((abs(lastTriangularizationIteration-iteration)) < minIterationsBeforeRetriangulation+min([5 round(iteration/50)]))
        reTriangulationIsNeeded = false;
    end
        
    % Always triangulate if this is the first iteration
    if (iteration == 1)
       reTriangulationIsNeeded = true;
       triangularizationTriggerEvent = '1st iteration';
    end
        
end

