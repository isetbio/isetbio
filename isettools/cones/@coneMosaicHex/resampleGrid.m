function resampleGrid(obj, resamplingFactor)
%  Sample original rectmosaic using hex grid sampled at resamplingFactor
%
% Syntax:
%   resampleGrid(obj, resamplingFactor)
%
% Description:
%    Sample the original rectangular mosaic using a hex grid sampled at the
%    passed resamplingFactor.
%
% Inputs:
%    obj              - The cone mosaic hex object
%    resamplingFactor - The resampling factor
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/21/18  jnm  Formatting

    % Restore original state
    obj.restoreOriginalResState();
    
    % Compute hex grid nodes
    obj.coneLocsHexGrid = computeHexGridNodes(obj);
    
    % Decrease patternSampleSize and increase the mosaicSize both by the
    % resamplingFactor
    obj.resamplingFactor = resamplingFactor;
    obj.patternSampleSize = obj.patternSampleSize / obj.resamplingFactor;
    obj.mosaicSize = obj.mosaicSize * obj.resamplingFactor;
    
    % Sample the hex grid at the nodes of the high res rectangular grid
    obj.pattern = rectSampledHexPattern(obj);
end

function hexLocs = computeHexGridNodes(obj)
% Compute minimum cone spacing (in microns)
%
% Syntax:
%   hexLocs = computeHexGridNodes(obj)
%
% Description:
%    Compute the minimum cone spacing of the provided cone mosaic hex, with
%    microns for the units.
%
% Inputs:
%    obj     - Cone mosaic hex object
%
% Outputs:
%    hexLocs - The minimum cone spacing, in microns.
%
% Optional key/value pairs:
%    None.
%
    obj.lambdaMin = minConeSpacing(obj);
    obj.lambdaMid = midConeSpacing(obj);
    grid.lambdaMin = obj.lambdaMin;
    grid.lambdaMid = obj.lambdaMid;
    if (~isempty(obj.customLambda))
        grid.lambdaMid = obj.customLambda;
        grid.lambdaMin = obj.customLambda;
    end

    grid.coneSpacingFunction = @coneSpacingFunction;
    grid.domainFunction = @ellipticalDomainFunction;
    grid.center = obj.center * 1e6;
    grid.rotationAngle = obj.rotationDegs / 180 * pi;
    grid.width = obj.width * 1e6;
    grid.height = obj.height * 1e6;
    grid.radius = obj.marginF * sqrt(2) * ...
        max([grid.width / 2, grid.height / 2]);
    grid.ellipseAxes = determineEllipseAxesLength(grid.radius);
    grid.borderTolerance = 0.001 * obj.lambdaMin;
    
    if (obj.eccBasedConeDensity)
        hexLocs = generateConePositionsOnVaryingDensityGrid(obj, grid);
    else
        grid.domainFunction = @circularDomainFunction;
        hexLocs = generateConePositionsOnConstantDensityGrid(grid);
    end

    % make sure we that the most foveal cone is at (0, 0)
    tmpHexLocs = bsxfun(@minus, hexLocs, obj.center);
    [~, fovealConeIndex] = min(sum(tmpHexLocs .^ 2, 2));
    xy0 = squeeze(tmpHexLocs(fovealConeIndex, :));
    tmpHexLocs = bsxfun(@minus, tmpHexLocs, xy0);
    hexLocs = bsxfun(@plus, tmpHexLocs, obj.center);

    % The cones within the rect mosaic extent
    mosaicRangeX = grid.center(1) + grid.width / 2 * [-1 1];
    mosaicRangeY = grid.center(2) + grid.height / 2 * [-1 1];
    if (obj.eccBasedConeDensity)
        indices = find( ...
            (hexLocs(:, 1) >= mosaicRangeX(1) + obj.lambdaMin / 4 + ...
            eps(mosaicRangeX(1) + obj.lambdaMin / 4)) & ...
            (hexLocs(:, 1) <= mosaicRangeX(2) - obj.lambdaMin / 4 - ...
            eps(mosaicRangeX(2) - obj.lambdaMin / 4)) & ... 
            (hexLocs(:, 2) >= mosaicRangeY(1) + obj.lambdaMin / 4 + ...
            eps(mosaicRangeY(1) + obj.lambdaMin / 4)) & ...
            (hexLocs(:, 2) <= mosaicRangeY(2) - obj.lambdaMin / 4 - ...
            eps(mosaicRangeY(2) - obj.lambdaMin / 4)) );
    else
        radii = sqrt(sum(hexLocs .^ 2, 2));
        indices = find( ...
            (hexLocs(:, 1) >= mosaicRangeX(1) + obj.lambdaMin / 4 + ...
            eps(mosaicRangeX(1) + obj.lambdaMin / 4)) & ...
            (hexLocs(:, 1) <= mosaicRangeX(2) - obj.lambdaMin / 4 - ...
            eps(mosaicRangeX(2) - obj.lambdaMin / 4)) & ... 
            (hexLocs(:, 2) >= mosaicRangeY(1) + obj.lambdaMin / 4 + ...
            eps(mosaicRangeX(1) + obj.lambdaMin / 4)) & ...
            (hexLocs(:, 2) <= mosaicRangeY(2) - obj.lambdaMin / 4 - ...
            eps(mosaicRangeX(2) - obj.lambdaMin / 4)) & ...
            (radii <= grid.radius - eps(grid.radius)) );
    end
    hexLocs = hexLocs(indices, :);

    % Return positions in meters
    hexLocs = hexLocs * 1e-6;
end

function conePositions = generateConePositionsOnVaryingDensityGrid(obj, ...
    gridParams) 
% Generate the cone positions on a grid of varying density
%
% Syntax:
%   conePositions = generateConePositionsOnVaryingDensityGrid(obj, ...
%       gridParams)
%
% Description:
%    Generate the cone positions on a grid of varying density
%
% Inputs:
%    obj           - The cone mosaic hex object
%    gridParams    - Struct containing grid parameters
%
% Outputs:
%    conePositions - The calculated positions of the cones.
%
% Optional key/value pairs:
%    None.
%
    % First generate perfect grid
    conePositions = generateConePositionsOnPerfectGrid(...
        gridParams.center, gridParams.radius, gridParams.lambdaMin, ...
        gridParams.rotationAngle);

    % Remove cones outside the desired region by applying the provided
    % domain function
    d = feval(gridParams.domainFunction, conePositions, ...
        gridParams.center, gridParams.radius, gridParams.ellipseAxes);
    conePositions = conePositions(d < gridParams.borderTolerance, :);
    
    if (obj.saveLatticeAdjustmentProgression)
        obj.initialLattice = conePositions * 1e-6;
        obj.latticeAdjustmentSteps = [];
    end
    
    % sample probabilistically according to coneSpacingFunction
    coneSeparations = feval(gridParams.coneSpacingFunction, conePositions);
    normalizedConeSeparations = coneSeparations / gridParams.lambdaMin;
    fudgeFactor = 1.0;
    densityP = fudgeFactor * (1 ./ normalizedConeSeparations) .^ 2;
    
    % Remove cones accordingly
    fixedConePositionsRadiusInCones = 1;
    radii = sqrt(sum(conePositions.^2,2));
    
    keptConeIndices = find(...
        (rand(size(conePositions, 1), 1) < densityP) | ...
        ((radii < fixedConePositionsRadiusInCones*gridParams.lambdaMin)) );
 
    conePositions = conePositions(keptConeIndices, :);

    % Add jitter
    beginWithJitteredPositions = false;
    if (beginWithJitteredPositions)
        conePositions = conePositions + randn(size(conePositions)) * ...
            gridParams.lambdaMin / 6;
    end
    
    if (obj.saveLatticeAdjustmentProgression)
        obj.latticeAdjustmentSteps(1, :, :) = conePositions * 1e-6;
    end
    % Iteratively adjust the grid for a smooth coverage of the space
    conePositions = smoothGrid(obj, conePositions, gridParams);
end

function conePositions = smoothGrid(obj, conePositions, gridParams)
% Iteratively adjust the grid for a smooth coverage of the space
%
% Syntax:
%   conePositions = smoothGrid(obj, conePositions, gridParams)
%
% Description:
%    Iteratively adjust the grid to ensure a smooth coverage of the space.
%
% Inputs:
%    obj           - The cone mosaic hex object
%    conePositions - The position of the cones
%    gridParams    - A struct containing the grid parameters
%
% Outputs:
%    conePositions - The modified cone positions.
%
% Optional key/value pairs:
%    None.
%
    % Convergence parameters

%     s = struct(...
%         'positionalToleranceF', obj.latticeAdjustmentPositionalToleranceF, ...
%         'DelaunayToleranceF', obj.latticeAdjustmentDelaunayToleranceF, ...
%         'maxGridAdjustmentIterations', obj.maxGridAdjustmentIterations, ...
%         'gridParams', gridParams, ...
%         'conePositions', conePositions);
%     save('s.mat', 's');
%     pwd
%     pause
    
    positionalDiffTolerance = obj.latticeAdjustmentPositionalToleranceF ...
        * gridParams.lambdaMin;
    deps = sqrt(eps) * gridParams.lambdaMin;

    deltaT = 0.2;
    dTolerance = obj.latticeAdjustmentDelaunayToleranceF * ...
        gridParams.lambdaMin;

    % Initialize convergence
    oldConePositions = inf;
    forceMagnitudes = [];

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');

    % Number of cones to be adjusted
    conesNum = size(conePositions, 1);

    %fixedConePositionsRadiusInCones = 10;
    %fixedConesPositions = find(sqrt(sum(conePositions.^2,2))< fixedConePositionsRadiusInCones*gridParams.lambdaMin);
    %fprintf('Fixing positions of %d (out of %d)central cones\n', numel(fixedConesPositions), size(conePositions,1))
    fixedConesPositions = [];
    
    % Iteratively adjust the cone positions until the forces between nodes
    % (conePositions) reach equlibrium.
    notConverged = true;
    terminateAdjustment = 0;
    nextQueryIteration = obj.queryGridAdjustmentIterations;
    iteration = 0;
    triangulationIndex = 0;
    triangulationIndex2 = 0;
    
    thresholdFractionalPositionChange = 0.9;
    
    tic
    while (notConverged) && (iteration <= obj.maxGridAdjustmentIterations) && (terminateAdjustment == 0)
        iteration = iteration + 1;

        if (obj.maxGridAdjustmentIterations < 100)
            fprintf('\nHex grid adjustment: on iteration %d ... ', ...
                iteration-1);
        else
            if (mod(iteration,50) == 1)
                fprintf('\nHex grid adjustment: on iteration %d ...', ...
                    iteration);
            end
        end


        % compute cone positional diffs
        positionalDiffs = sqrt(...
            sum((conePositions-oldConePositions) .^ 2, 2));

        doDelaunaynTriangulation = false;
        if (max(positionalDiffs) > positionalDiffTolerance)
            doDelaunaynTriangulation = true;
        else
            % check if we need to update the mosaic triangulation (determine springs between cone neigbors)
            coneSeparations = feval(gridParams.coneSpacingFunction, conePositions);
            
            if ((max(positionalDiffs ./ coneSeparations) > thresholdFractionalPositionChange))
                doDelaunaynTriangulation = true;
                triangulationIndex2 = triangulationIndex2+1;
                figure(99); clf;
                plot(coneSeparations, positionalDiffs ./ coneSeparations, 'k.');
                hold on;
                plot([1 7], thresholdFractionalPositionChange*[1 1], 'r--');
                hold off;
                title(sprintf('iteration: %d', iteration));
                xlabel('cone separation');
                ylabel('positionalDiff/coneSeparation');
                drawnow
            end
        end
        
        if (doDelaunaynTriangulation)
            triangulationIndex = triangulationIndex + 1;
            
            % save old come positions
            oldConePositions = conePositions;

            % Perform new Delaunay triangulation to determine the updated
            % topology of the truss. To save computing time, we
            % re-triangulate only when we exceed the
            % positionalDiffTolerance
            triangleConeIndices = delaunayn(conePositions);

            % Compute the centroids of all triangles
            centroidPositions = (conePositions(...
                triangleConeIndices(:, 1), :) + conePositions(...
                triangleConeIndices(:, 2), :) + conePositions(...
                triangleConeIndices(:, 3), :)) / 3;

            % Remove centroids outside the desired region by applying the
            % signed distance function
            d = feval(gridParams.domainFunction, centroidPositions, ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes);
            triangleConeIndices = triangleConeIndices(d < ...
                gridParams.borderTolerance, :);

            % Create a list of the unique springs (each spring connecting 2
            % cones)
            % triangleConeIndices is an [M x 3] matrix the m-th row
            % contains indices to the 3 cones that define the triangle
            springs = [triangleConeIndices(:, [1, 2]);
                triangleConeIndices(:, [1, 3]);
                triangleConeIndices(:, [2, 3])];
            springs = unique(sort(springs, 2), 'rows');

            % find all springs connected to this cone
            for coneIndex = 1:conesNum
                springIndices{coneIndex} = find(...
                    (springs(:, 1) == coneIndex) | ...
                    (springs(:, 2) == coneIndex));
            end
        end % if (doDelaunaynTriangulation)
        
        % Compute spring vectors
        springVectors =  conePositions(springs(:, 1), :) - ...
            conePositions(springs(:, 2), :);
        % their centers
        springCenters = (conePositions(springs(:, 1), :) + ...
            conePositions(springs(:, 2), :)) / 2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2, 2));
        
        % Compute desired spring lengths. This is done by evaluating the
        % passed coneDistance function at the spring centers.
        desiredSpringLengths = feval(gridParams.coneSpacingFunction, ...
            springCenters);

        % Normalize spring lengths
        normalizingFactor = sqrt(sum(springLengths .^ 2) / ...
            sum(desiredSpringLengths .^ 2));
        desiredSpringLengths = desiredSpringLengths * normalizingFactor;

        % Compute spring forces, Force(springLengths, desiredSpringLengths)
        % Force(springLengths, desiredSpringLengths) should be positive
        % when springLengths is near the desiredSpringLengths, which can be
        % achieved by choosing desiredSpringLengths slightly larger than
        % the length we actually desire. Here, we set this to be 1.2
        gain = 1.2;
        springForces = max(gain * desiredSpringLengths - springLengths, 0);

        % compute x, y-components of forces on each of the springs
        springForceXYcomponents = abs(springForces ./ springLengths * ...
            [1, 1] .* springVectors);

        % Compute net forces on each cone
        netForceVectors = zeros(conesNum, 2);
        parfor coneIndex = 1:conesNum
            % compute net force from all connected springs
            deltaPos = -bsxfun(@minus, springCenters(...
                springIndices{coneIndex}, :), conePositions(coneIndex, :));
            netForceVectors(coneIndex, :) = sum(sign(deltaPos) .* ...
                springForceXYcomponents(springIndices{coneIndex}, :), 1);
        end

        % force at all fixed cone positions must be 0
        netForceVectors(1:size(fixedConesPositions, 1), :) = 0;
        
        % Save force magnitudes
        % forceMagnitudes(iteration, :) = ...
        %    sqrt(sum(netForceVectors .^ 2, 2)) / gridParams.lambdaMin;
        
        % update cone positions according to netForceVectors
        conePositions = conePositions + deltaT * netForceVectors;
        
        % Find any points that lie outside the domain boundary
        d = feval(gridParams.domainFunction, conePositions, ...
            gridParams.center, gridParams.radius, gridParams.ellipseAxes);
        outsideBoundaryIndices = d > 0;
        
        % And project them back to the domain
        if (~isempty(outsideBoundaryIndices))
            % Compute numerical gradient along x-positions
            dXgradient = (feval(gridParams.domainFunction, ...
                [conePositions(outsideBoundaryIndices, 1) + deps, ...
                conePositions(outsideBoundaryIndices, 2)], ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                deps;
            dYgradient = (feval(gridParams.domainFunction, ...
                [conePositions(outsideBoundaryIndices, 1), ...
                conePositions(outsideBoundaryIndices, 2)+deps], ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                deps;

            % Project these points back to boundary
            conePositions(outsideBoundaryIndices, :) = ...
                conePositions(outsideBoundaryIndices, :) - ...
                [d(outsideBoundaryIndices) .* dXgradient, ...
                d(outsideBoundaryIndices) .* dYgradient];
        end

        % Check if all interior nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(...
            d < -gridParams.borderTolerance, :) .^2 , 2));
        if max(movementAmplitudes) < dTolerance, notConverged = false; end

        if (obj.saveLatticeAdjustmentProgression)
            obj.latticeAdjustmentSteps(...
                size(obj.latticeAdjustmentSteps, 1) + 1, :, :) = ...
                conePositions * 1e-6;
        end 
        
        % check whether we need to ask user whether to continue or not
        if (iteration == nextQueryIteration) % (mod(iteration,obj.queryGridAdjustmentIterations) == 0)
            visualizeLatticeState(obj, conePositions, iteration);
            hoursLapsed = toc/60/60;
            qString = sprintf('\n[at iter %d after %2.2f hours] Terminate adjusting (1) or continue (0)', iteration, hoursLapsed);
            terminateAdjustment = queryUserWithDefault(qString, 0);
            if (terminateAdjustment == 0)
                possibleIterationIntervals = obj.queryGridAdjustmentIterations * [0.5 1 2 5 10];
                nextQueryIterationString = sprintf('Next query iterations. 1=%2.0f 2=%2.0f 3=%2.0f 4=%2.0f 5=%2.0f : ', ...
                    possibleIterationIntervals(1), possibleIterationIntervals(2), possibleIterationIntervals(3), possibleIterationIntervals(4), possibleIterationIntervals(5));
                nextQueryIteration = queryUserWithDefault(nextQueryIterationString, 1);
                if (nextQueryIteration < 1)
                    nextQueryIteration = 1;
                elseif (nextQueryIteration > numel(possibleIterationIntervals))
                    nextQueryIteration = numel(possibleIterationIntervals);
                end
                nextQueryIteration = iteration+possibleIterationIntervals(nextQueryIteration);
                fprintf('Will ask again at iteration %d.\n', nextQueryIteration);
                tic
            else
                fprintf('Terminating adjustment at user request\n');
            end
        else
            if (~isinf(obj.maxGridAdjustmentIterations)) && (mod(iteration,obj.visualizationUpdateIterations) == 0)
                visualizeLatticeState(obj, conePositions, iteration);
            end
        end
        
    end % while (notConverged) && (iteration < obj.maxGridAdjustmentIterations)
    
    fprintf('\nHex grid smoothing finished.');
    if (iteration > obj.maxGridAdjustmentIterations) 
        fprintf('\nDid not converge, but exceeded max number of iterations (%d).', ...
            obj.maxGridAdjustmentIterations);
        fprintf('\nMax(movement) in last iteration: %2.6f, Tolerange: %2.6f\n', ...
            max(movementAmplitudes), obj.latticeAdjustmentDelaunayToleranceF);
    else
        fprintf('Converged after %d iterations.\n', iteration);
    end
    fprintf('Number of Delayun triangulations: %d\n', triangulationIndex);
    fprintf('Number of Delayun triangulations due to fractional movement: %d\n', triangulationIndex2);
    
    % Turn back on Delaunay triangularization warning
    warning('on', 'MATLAB:qhullmx:InternalWarning');
end

function conePositions = generateConePositionsOnConstantDensityGrid(...
    gridParams)
% Generate the cone positions on a grid with a constant density
%
% Syntax:
%   conePositions = generateConePositionsOnConstantDensityGrid(gridParams)
%
% Description:
%    Generate the positions of the cones on a grid that has a constant
%    density, using the provided parameters.
%
% Inputs:
%    gridParams    - A struct containing the grid params
%
% Outputs:
%    conePositions - The calculated cone positions.
%
% Optional key/value pairs:
%    None.
%
    conePositions = generateConePositionsOnPerfectGrid(...
        gridParams.center, gridParams.radius, gridParams.lambdaMid, ...
        gridParams.rotationAngle);
    % Remove cones outside of the desired region by applying the passed
    % domain function
    d = feval(gridParams.domainFunction, conePositions, ...
        gridParams.center, gridParams.radius, gridParams.ellipseAxes);
    conePositions = conePositions(d < gridParams.borderTolerance, :);
end

function conePositions = generateConePositionsOnPerfectGrid(...
    center, radius, lambda, rotationAngle)
% Generate the cone positions on a perfect grid
%
% Syntax:
%   conePositions = generateConePositionsOnPerfectGrid(center, radius, ...
%   lambda, rotationAngle)
%
% Description:
%    Generate the position of the cones on a perfect grid.
%
% Inputs:
%    center        - The (X, Y) location of the center of the grid
%    radius        - The radius of the grid (1/2 number of rows & cols)
%    lambda        - The distance between the cones
%    rotationAngle - The angle of rotation for the grid
%
% Outputs:
%    conePositions - The calculated cone positions
%
% Optional key/value pairs:
%    None.
%
    rows = 2 * radius;
    cols = 2 * radius;
    conePositions = computeHexGrid(rows, cols, lambda, rotationAngle);
    conePositions(:, 1) = conePositions(:, 1) + center(1);
    conePositions(:, 2) = conePositions(:, 2) + center(2);
end

function lambda = midConeSpacing(obj)
% Calculate the spacing between cones, in microns
%
% Syntax:
%   lambda = midConeSpacing(obj)
%
% Description:
%    Calculate the spacing between cones on a grid.
%
% Inputs:
%    obj    - The cone mosaic hex object
%
% Outputs:
%    lambda - The calculated spacing between cones in microns
%
% Optional key/value pairs:
%    None.
%
    midX = obj.center(1);
    midY = obj.center(2);
    eccentricityInMeters = sqrt(midX ^ 2 + midY ^ 2);
    ang = atan2(midY, midX) / pi * 180;
    [coneSpacingInMeters, aperture, density] = coneSizeReadData(...
        'eccentricity', eccentricityInMeters, 'angle', ang);
    lambda = coneSpacingInMeters * 1e6;  % in microns
end

function lambda = minConeSpacing(obj)
% Calculate the minimum space between cones, in microns
%
% Syntax:
%   lambda = minConeSpacing(obj)
%
% Description:
%    Calculate the minimum spacing between the cones, in microns.
%
% Inputs:
%    obj    - The cone mosaic hex object
%
% Outputs:
%    lambda - The calculated spacing between cones in microns
%
% Optional key/value pairs:
%    None.
%
    mosaicRangeX = obj.center(1) + obj.width / 2 * [-1 1];
    mosaicRangeY = obj.center(2) + obj.height / 2 * [-1 1];
    minX = min(abs(mosaicRangeX));
    if (prod(mosaicRangeX) < 0), minX = 0; end
    minY = min(abs(mosaicRangeY));
    if (prod(mosaicRangeY) < 0), minY = 0; end
    eccentricityInMeters = sqrt(minX ^ 2 + minY ^ 2);
    ang = atan2(minY, minX) / pi * 180;
    [coneSpacingInMeters, aperture, density] = coneSizeReadData(...
        'eccentricity', eccentricityInMeters, 'angle', ang);
    lambda = coneSpacingInMeters * 1e6;  % in microns
end

function hexLocs = computeHexGrid(rows, cols, lambda, rotationAngle)
% Compute the hex grid
%
% Syntax:
%   hexLocs = computeHexGrid(rows, cols, lambda, rotationAngle)
%
% Description:
%    Compute the hex grid.
%
% Inputs:
%    rows          - The number of rows
%    cols          - The number of columns
%    lambda        - The space between, to help calculate density
%    rotationAngle - The angle at which to rotate the grid
%
% Outputs:
%    hexLocs       - The location of the hexes
%
% Optional key/value pairs:
%    None.
%
    scaleF = sqrt(3) / 2;
    extraCols = round(cols / scaleF) - cols;
    rectXaxis2 = (1:(cols + extraCols));
    [X2, Y2] = meshgrid(rectXaxis2, 1:rows);

    X2 = X2 * scaleF ;
    for iCol = 1:size(Y2, 2)
        Y2(:, iCol) = Y2(:, iCol) - mod(iCol - 1, 2) * 0.5;
    end

    % Scale to get correct density
    X2 = X2 * lambda;
    Y2 = Y2 * lambda;
    marginInConePositions = 0.1;
    indicesToKeep = (X2 >= -marginInConePositions) & ...
                    (X2 <= cols+marginInConePositions) &...
                    (Y2 >= -marginInConePositions) & ...
                    (Y2 <= rows+marginInConePositions);
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
    hexLocs = [xHex(:) - mean(xHex(:)) yHex(:) - mean(yHex(:))];

    % rotate
    R = [cos(rotationAngle) -sin(rotationAngle);
         sin(rotationAngle) cos(rotationAngle)];
    hexLocs = (R * (hexLocs)')';
end

function distances = circularDomainFunction(conePositions, center, ...
    radius, nullVar)
% Calculate distances for a ciruclar domain function
%
% Syntax:
%   distances = circularDomainFunction(conePositions, center, radius)
%
% Description:
%    Calculate the circular domain function distances.
%
% Inputs:
%    conePositions - The cone grid positions
%    center        - The center of the grid
%    radius        - The radius of the grid
%
% Outputs:
%    distances     - Distances in a circular domain
%
% Optional key/value pairs:
%    None.
%
    %  points with positive distance will be excluded
    radii = sqrt((conePositions(:, 1) - center(1)) .^ 2 + ...
        (conePositions(:, 2) - center(2)) .^ 2);
    distances = -(radius - radii);
end

function distances = ellipticalDomainFunction(conePositions, center, ...
    radius, ellipseAxes)
% Calculate distances for an elliptical domain function
%
% Syntax:
%   distances = ellipticalDomainFunction(conePositions, center, radius, ...
%       ellipseAxes)
%
% Description:
%    Calcularte the elliptical domain function distances.
%
% Inputs:
%    conePositions - The cone grid positions
%    center        - The center of the grid
%    radius        - The radius of the grid
%    ellipseAxes   - The ellipse's axes
%
% Outputs:
%    distances     - The distances in the elliptical domain
%
% Optional key/value pairs:
%    None.
%
    %  points with positive distance will be excluded
    xx = conePositions(:, 1) - center(1);
    yy = conePositions(:, 2) - center(2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = -(radius - radii);
end

function ellipseAxes = determineEllipseAxesLength(mosaicHalfFOVmicrons)
% compute ellipse axes lengths
%
% Syntax:
%    ellipseAxes = determineEllipseAxesLength(mosaicHalfFOVmicrons)
%
% Description:
%    Calculate the length of the ellipse's axes.
%
% Inputs:
%    mosaicHalfFOVmicrons - Half o the mosaic's field of view, in microns
%
% Outputs:
%    ellipseAxes          - The length of the ellipse's axes
%
% Optional key/value pairs:
%    None.
%
    largestXYspacing = coneSizeReadData('eccentricity', ...
        1e-6 * mosaicHalfFOVmicrons * [1 1], 'angle', [0 90]);
    if (largestXYspacing(1) < largestXYspacing(2))
        xa = 1;
        ya = largestXYspacing(2) / largestXYspacing(1);
    else
        xa = largestXYspacing(1) / largestXYspacing(2);
        ya = 1;
    end
    ellipseAxes = [xa ya] .^ 2;
end

function [coneSpacingInMicrons, eccentricitiesInMicrons] = ...
    coneSpacingFunction(conePositions)
% Calculate the cone spacing
%
% Syntax:
%   [coneSpacingInMicrons, eccentricitiesInMicrons] = ...
%       coneSpacingFunction(conePositions)
%
% Description:
%    Calculate the cone spacing and eccentricities, both in microns.
%
% Inputs:
%    conePositions           - The position of the cones
%
% Outputs:
%    coneSpacingInMicrons    - The spacing between cones, in microns.
%    eccentricitiesInMicrons - The cone eccentricities, in microns.
%
% Optional key/value pairs:
%    None.
%
    eccentricitiesInMicrons = sqrt(sum(conePositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(conePositions(:, 2), conePositions(:, 1)) / pi * 180;
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6;
end

function pattern = rectSampledHexPattern(obj)

% Create a high resolution rectangular mosaic pattern
%
% Syntax:
%   pattern = rectSampledHexPattern(obj)
%
% Description:
%    create a high resolution rectangular mosaic pattern using the provided
%    cone mosaic hex object.
%
% Inputs:
%    obj     - A cone mosaic hex object
%
% Outputs:
%    pattern - The created rectangular mosaic pattern 
%
% Optional key/value pairs:
%    None.
%
    fprintf('\nResampling grid. Please wait ... ');
    % High resolution grid
    xRectHiRes = (1:obj.cols) * obj.patternSampleSize(1);
    xRectHiRes = xRectHiRes - mean(xRectHiRes);
    yRectHiRes = (1:obj.rows) * obj.patternSampleSize(2);
    yRectHiRes = yRectHiRes - mean(yRectHiRes);
    [xx, yy] = meshgrid(xRectHiRes, yRectHiRes);

    % Match spatial extent
    xRectOriginal = (1:size(obj.patternOriginatingRectGrid, 2)) * ...
        obj.patternSampleSizeOriginatingRectGrid(1);
    xRectOriginal = xRectOriginal - mean(xRectOriginal);
    yRectOriginal = (1:size(obj.patternOriginatingRectGrid, 1)) * ...
        obj.patternSampleSizeOriginatingRectGrid(2);
    yRectOriginal = yRectOriginal - mean(yRectOriginal);
    xRectOriginal = xRectOriginal / max(xRectOriginal) * max(xRectHiRes);
    yRectOriginal = yRectOriginal / max(yRectOriginal) * max(yRectHiRes);
    [xxx, yyy] = meshgrid(xRectOriginal, yRectOriginal);

    % Generate the high res mosaic pattern
    pattern = zeros(numel(yRectHiRes), numel(xRectHiRes)) + 1;

    % Determine the closest cone type in the originating  grid
    [~, I] = pdist2([xxx(:) yyy(:)], bsxfun(@minus, ...
        obj.coneLocsHexGrid, obj.center), 'euclidean', 'Smallest', 1);

    % Determine the closest cone location in the high-res grid
    [~, II] = pdist2([xx(:) yy(:)], bsxfun(@minus, ...
        obj.coneLocsHexGrid, obj.center), 'euclidean', 'Smallest', 1);

    % That's our cone!
    pattern(ind2sub(size(obj.pattern), II)) = ...
        obj.patternOriginatingRectGrid(I);

    % Make cones outside of the desired FOV null
    outsideBoundaryIndices = find((abs(xx(:) - obj.center(1)) > ...
        obj.width / 2) & (abs(yy(:) - obj.center(2)) > obj.height / 2) );
    pattern(ind2sub(size(obj.pattern), outsideBoundaryIndices)) = 0;

    fprintf('Done !\n');

    % Show patterns (only for debugging purposes)
    debugPlots = false;
    if (debugPlots)
        cmap = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
        figure(10);
        clf;
        imagesc((obj.patternOriginatingRectGrid - 1) / 3);
        axis 'equal';
        axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'

        figure(11);
        clf;
        imagesc((pattern - 1) / 3);
        axis 'equal';
        axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'
    end
end

function visualizeLatticeState(obj, conePositions, iteration)
    qDist = computeQuality(conePositions);
    threshold = 150;
    r = sqrt(sum(conePositions.^2,2));
    idx = find(r < threshold);
    conePositions = conePositions(idx,:);
    hFig = figure(111); clf;
    set(hFig,'Position', [10 10 1650 950]);
    subplot('Position', [0.01 0.02 0.68 0.97]);
    plot(conePositions(:,1), conePositions(:,2), 'ko', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 4);
    set(gca, 'XLim', threshold*[-1 1], 'YLim', threshold*[-1 1], 'XTick', [], 'YTick', []);
    axis 'square';
    subplot('Position', [0.71 0.04 0.28 0.94]);
    plotQuality(qDist);
    title(sprintf('iteration %d', iteration))
end

function plotQuality(qDist)
    qLims = [0.5 1.005]; qBins = [0.3:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 max(counts)], 'XTick', [0.1:0.1:1.0],  'FontSize', 16);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('count', 'FontSize', 16);
end

function q = computeQuality(coneLocs)
    
    triangles = delaunay(squeeze(coneLocs(:,1)),squeeze(coneLocs(:,2)));
    
    trianglesNum = size(triangles,1);
    X = coneLocs(:,1);
    Y = coneLocs(:,2);
    
    q = zeros(1,trianglesNum);
    for triangleIndex = 1:trianglesNum
        for node = 1:3
            x(node) = X(triangles(triangleIndex,node));
            y(node) = Y(triangles(triangleIndex,node));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        q(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end
end

function inputVal = queryUserWithDefault(prompt,defaultVal)
    if (ischar(defaultVal))
        inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
    else
        inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
    end
    if (isempty(inputVal))
        inputVal = defaultVal;
    end
end