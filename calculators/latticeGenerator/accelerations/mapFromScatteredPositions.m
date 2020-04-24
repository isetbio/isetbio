function [densityMap, support] = mapFromScatteredPositions(rfPositions, values, sampling)

    % Now interpolate to a regular 2D grid
    F = scatteredInterpolant(rfPositions(:,1),rfPositions(:,2), values(:), 'linear');
    
    % Compute the support for the density map
    if (strcmp(sampling.scale, 'log'))
        if (sampling.minPos == 0)
            error('minPos cannot be 0 with a log spacing');
        end
        xAxis = logspace(...
            log10(sampling.minPos), ...
            log10(sampling.maxPos), ...
            sampling.intervals);
        xAxis = [-fliplr(xAxis(2:end)) 0 xAxis];
    else
        xAxis = linspace(...
            sampling.minPos, ...
            sampling.maxPos, ...
            sampling.intervals);
        xAxis = [-fliplr(xAxis(2:end)) xAxis];
    end
    
    [X,Y] = meshgrid(xAxis, xAxis);

    % Evaluate the interpolant function at the support
    densityMap = F(X, Y);
    support(:,:,1) = X;
    support(:,:,2) = Y;
end

