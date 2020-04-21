function [densityMap, support] = densityMapFromPositions(rfPositions, densitySampling)
    neighborsNum = 5;

    p = pdist2(rfPositions, rfPositions, 'euclidean', 'Smallest', neighborsNum+1);
    % exclude the points themselves which have 0 distance
    p = p(2:end,:);
    
    % Compute the median spacing among the neibhors
    spacings = median(p,1);
    densities = WatsonRGCModel.densityFromSpacing(spacings);
    
    % Now interpolate to a regular 2D grid
    F = scatteredInterpolant(rfPositions(:,1),rfPositions(:,2), densities(:), 'linear');
    
    % Compute the support for the density map
    if (strcmp(densitySampling.scale, 'log'))
        xAxis = logspace(...
            log10(densitySampling.minPos), ...
            log10(densitySampling.maxPos), ...
            densitySampling.intervals);
    else
        xAxis = linspace(...
            densitySampling.minPos, ...
            densitySampling.maxPos, ...
            densitySampling.intervals);
    end
    xAxis = [-fliplr(xAxis(2:end)) xAxis];
    [X,Y] = meshgrid(xAxis, xAxis);

    % Evaluate the interpolant function at the support
    densityMap = F(X, Y);
    support(:,:,1) = X;
    support(:,:,2) = Y;
end

