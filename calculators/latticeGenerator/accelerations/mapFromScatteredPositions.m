function [densityMap, support] = mapFromScatteredPositions(rfPositions, values, sampling)

    % Now interpolate to a regular 2D grid
    F = scatteredInterpolant(rfPositions(:,1),rfPositions(:,2), values(:), 'linear');
    
    % Compute the support for the density map
    if (strcmp(sampling.scale, 'log'))
        xAxis = logspace(...
            log10(sampling.minPos), ...
            log10(sampling.maxPos), ...
            sampling.intervals);
    else
        xAxis = linspace(...
            sampling.minPos, ...
            sampling.maxPos, ...
            sampling.intervals);
    end
    xAxis = [-fliplr(xAxis(2:end)) xAxis];
    [X,Y] = meshgrid(xAxis, xAxis);

    % Evaluate the interpolant function at the support
    densityMap = F(X, Y);
    support(:,:,1) = X;
    support(:,:,2) = Y;
end

