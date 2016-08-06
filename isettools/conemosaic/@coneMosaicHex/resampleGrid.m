function resampleGrid(obj, resamplingFactor)
%  Sample the original rectangular mosaic using a hex grid sampled at the passed resamplingFactor
%  
% NPC, ISETBIO TEAM, 2015

    % Restore original state
    obj.restoreOriginalResState();
    
    % Compute hex grid nodes
    obj.coneLocsHexGrid = computeHexGridNodes(obj);
    %obj.coneLocsHexGrid = computeHexGridNodesOriginal(obj);
    
    % Decrease patternSampleSize and increase mosaicSize both by the resamplingFactor
    obj.resamplingFactor = resamplingFactor;
    obj.patternSampleSize = obj.patternSampleSize / obj.resamplingFactor;
    obj.mosaicSize = obj.mosaicSize * obj.resamplingFactor;
    
    % Sample the hex grid at the nodes of the high res rectangular grid
    obj.pattern = rectSampledHexPattern(obj);
end

function hexLocs = computeHexGridNodes(obj)

    % Compute minimum cone spacing (in microns)
    obj.lambdaMin = minConeSpacing(obj);
    obj.lambdaMid = midConeSpacing(obj);
    grid.lambdaMin = obj.lambdaMin;
    grid.lambdaMid = obj.lambdaMid;
    grid.coneSpacingFunction = @coneSpacingFunction;
    grid.domainFunction = @circularDomainFunction;
    grid.center = obj.center*1e6;
    grid.width = obj.width*1e6;
    grid.height = obj.height*1e6;
    grid.radius = 1.1*sqrt((grid.width/2)^2 + (grid.height/2)^2);
    grid.borderTolerance = 0.001*obj.lambdaMin;
    
    if (obj.varyingDensity)
        hexLocs = generateConePositionsOnVaryingDensityGrid(grid);
    else
        hexLocs = generateConePositionsOnConstantDensityGrid(grid);
    end
    
    % The cones within the rect mosaic extent
    mosaicRangeX = grid.center(1) + grid.width/2*[-1 1];
    mosaicRangeY = grid.center(2) + grid.height/2*[-1 1];
    indices = find( ...
        hexLocs(:,1) >= mosaicRangeX(1)-obj.lambdaMin/4 & ...
        hexLocs(:,1) <= mosaicRangeX(2)+obj.lambdaMin/4 & ... 
        hexLocs(:,2) >= mosaicRangeY(1)-obj.lambdaMin/4 & ...
        hexLocs(:,2) <= mosaicRangeY(2)+obj.lambdaMin/4 );
    hexLocs = hexLocs(indices,:);
    
    % Return positions in meters
    hexLocs = hexLocs * 1e-6; 
end

function conePositions = generateConePositionsOnVaryingDensityGrid(gridParams)

    % generate mosaic of highest density
    conePositions = generateInitialConePositionsOnVaryingDensityGrid(gridParams);
end

function conePositions = generateInitialConePositionsOnVaryingDensityGrid(gridParams)
    
    % First generate perfect grid
    conePositions = generateConePositionsOnPerfectGrid(gridParams.center, gridParams.radius, gridParams.lambdaMin);
    
    % sample probabilistically according to coneSpacingFunction
    coneSeparations = feval(gridParams.coneSpacingFunction, conePositions);
    normalizedConeSeparations = coneSeparations/gridParams.lambdaMin;
    densityP = 0.956 * 2/sqrt(3) * (1 ./ normalizedConeSeparations).^2;
    
    showConeSeparations = false;
    if (showConeSeparations)
        figure(9); clf;
        scatter3(conePositions(:,1), conePositions(:,2), coneSeparations);
    end
    
    keptConeIndices = find(rand(size(conePositions,1), 1) < densityP);
    conePositions = conePositions(keptConeIndices,:);
    
    % Remove cones outside the desired region by applying the passed domainfunction
    d = feval(gridParams.domainFunction, conePositions, gridParams.center, gridParams.radius);
    conePositions = conePositions(d < gridParams.borderTolerance,:);
    
    % remove dublicate points
    [~, ia, ~] = unique(round(conePositions/gridParams.lambdaMin), 'rows');
    conePositions = conePositions(ia,:);

    % Add jitter
    beginWithJitteredPositions = false;
    if (beginWithJitteredPositions)
        conePositions = conePositions + randn(size(conePositions))*gridParams.lambdaMin/6;
    end
    
    % Iteratively adjust the grid for a smooth coverage of the space
    conePositions = smoothGrid(conePositions, gridParams);
end

% Iteratively adjust the grid for a smooth coverage of the space
function conePositions = smoothGrid(conePositions, gridParams)

    conesNum = size(conePositions,1);
%     continueCheck = input(sprintf('The mosaic will contain %d cones. Carry on with iterative spacing smoothing? [y/n]: ', conesNum), 's');
%     if (~strcmp(continueCheck, 'y'))
%         return;
%     end
    
    % Convergence parameters
    % fraction = 0.01;
    fraction = 0.5;
    positionalDiffTolerance = fraction * gridParams.lambdaMin;  
    deps = sqrt(eps)*gridParams.lambdaMin; 
    deltaT = 0.2;
    
    % fraction = 0.001;
    fraction = 0.05;
    dTolerance = fraction * gridParams.lambdaMin;
    
    % Initialize convergence
    oldConePositions = inf;
    forceMagnitudes = [];
    
    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    
    % Iteratively adjust the cone positions until the forces between nodes (conePositions) reach equlibrium.
    notConverged = true;
    iteration = 0;
    tic
    while (notConverged)
        iteration = iteration + 1;
        if (mod(iteration,50) == 1)
            fprintf('\nIteration: %d', iteration-1);
        end
        
        % compute cone positional diffs
        positionalDiffs = sqrt(sum((conePositions-oldConePositions).^2,2));
        
        % check if there are any large movements
        if (max(positionalDiffs) > positionalDiffTolerance)
            % save old come positions
            oldConePositions = conePositions;
            
            % Perform new Delaunay triangulation to determine the updated topology of the truss.
            % To save computing time, we re-triangulate only when we exceed the positionalDiffTolerance
            triangleConeIndices = delaunayn(conePositions);
            
            % Compute the centroids of all triangles
            centroidPositions = (conePositions(triangleConeIndices(:,1),:) + conePositions(triangleConeIndices(:,2),:) + conePositions(triangleConeIndices(:,3),:))/3; 
    
            % Remove centroids outside the desired region by applying the signed distance function
            d = feval(gridParams.domainFunction, centroidPositions, gridParams.center, gridParams.radius);
            triangleConeIndices = triangleConeIndices(d < gridParams.borderTolerance,:);
            
            % Create list of unique springs (each spring connecting 2 cones)
            % triangleConeIndices is an [M x 3] matrix the m-th row contains indices to the 3 cones that define the triangle
            springs = [triangleConeIndices(:, [1,2]); triangleConeIndices(:, [1,3]); triangleConeIndices(:, [2,3])];
            springs = unique(sort(springs, 2), 'rows');
            
            % find all springs connected to this cone
            for coneIndex = 1:conesNum
                springIndices{coneIndex} = find((springs(:,1) == coneIndex) | (springs(:,2) == coneIndex));
            end
        end % (max(positionalDiffs) > positionalDiffTolerance)
        
        % Compute spring vectors
        springVectors =  conePositions(springs(:,1),:) - conePositions(springs(:,2),:);
        % their centers
        springCenters = (conePositions(springs(:,1),:) + conePositions(springs(:,2),:))/2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2,2));
        
        % Compute desired spring lengths 
        % This is done by evaluating the passed coneDistance function at the spring centers
        desiredSpringLengths = feval(gridParams.coneSpacingFunction, springCenters);

        % Normalize spring lengths
        normalizingFactor = sqrt(sum(springLengths.^2)/sum(desiredSpringLengths.^2));
        desiredSpringLengths = desiredSpringLengths * normalizingFactor;
        
        % Compute spring forces, Force(springLengths, desiredSpringLengths)
        % Force(springLengths, desiredSpringLengths) should be positive when springLengths is near the desiredSpringLengths, 
        % which can be achieved by choosing desiredSpringLengths slightly larger than the length we actually desire 
        % Here, we set this to be 1.2
        gain = 1.2;
        springForces = max(gain * desiredSpringLengths - springLengths, 0);
        
        % compute x,y-components of forces on each of the springs
        springForceXYcomponents = abs(springForces ./ springLengths * [1,1] .* springVectors);
        
        % Compute net forces on each cone
        netForceVectors = zeros(conesNum,2);
        for coneIndex = 1:conesNum
            % compute net force from all connected springs
            deltaPos = -bsxfun(@minus, springCenters(springIndices{coneIndex},:), conePositions(coneIndex,:));
            netForceVectors(coneIndex,:) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{coneIndex},:), 1);   
        end
        
        % force at all fixed cone positions must be 0
        % netForceVectors(1:size(fixedConesPositions,1),:) = 0;
        forceMagnitudes(iteration,:) = sqrt(sum(netForceVectors.^2,2))/gridParams.lambdaMin;
        
        % update cone positions according to netForceVectors
        conePositions = conePositions + deltaT * netForceVectors;
        
        % Find any points that lie outside the domain boundary
        d = feval(gridParams.domainFunction, conePositions, gridParams.center, gridParams.radius);
        outsideBoundaryIndices = d > 0;
        
        % And project them back to the domain
        if (~isempty(outsideBoundaryIndices))
            % Compute numerical gradient along x-positions
            dXgradient = (feval(gridParams.domainFunction,[conePositions(outsideBoundaryIndices,1)+deps, conePositions(outsideBoundaryIndices,2)], gridParams.center, gridParams.radius) - d(outsideBoundaryIndices))/deps;
            dYgradient = (feval(gridParams.domainFunction,[conePositions(outsideBoundaryIndices,1), conePositions(outsideBoundaryIndices,2)+deps], gridParams.center, gridParams.radius) - d(outsideBoundaryIndices))/deps;

            % Project these points back to boundary
            conePositions(outsideBoundaryIndices,:) = conePositions(outsideBoundaryIndices,:) - ...
                [d(outsideBoundaryIndices).*dXgradient, d(outsideBoundaryIndices).*dYgradient]; 
        end
        
        % Check if all interior nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d < -gridParams.borderTolerance,:).^2,2));
        if max(movementAmplitudes) < dTolerance
            notConverged = false;
        end
        
    end % while (notConverged)
    fprintf('\nDone with iterative adjustment in %2.1f seconds\n', toc)
    
    % Turn back on Delaunay triangularization warning
    warning('on', 'MATLAB:qhullmx:InternalWarning');
end


function conePositions = generateConePositionsOnConstantDensityGrid(gridParams)
    conePositions = generateConePositionsOnPerfectGrid(gridParams.center, gridParams.radius, gridParams.lambdaMid);
    % Remove cones outside the desired region by applying the passed domainfunction
    d = feval(gridParams.domainFunction, conePositions, gridParams.center, gridParams.radius);
    conePositions = conePositions(d < gridParams.borderTolerance,:);
end

function conePositions = generateConePositionsOnPerfectGrid(center, radius, lambda)
    rows = 2*radius;
    cols = 2*radius;
    conePositions = computeHexGrid(rows, cols, lambda);
    conePositions(:,1) = conePositions(:,1) + center(1);
    conePositions(:,2) = conePositions(:,2) + center(2);
end

function lambda = midConeSpacing(obj)
    midX = obj.center(1);
    midY = obj.center(2);
    eccentricityInMeters = sqrt(midX^2 + midY^2);
    ang = atan2(midY, midX)/pi*180;
    [coneSpacingInMeters, aperture, density] = coneSize(eccentricityInMeters,ang);
    lambda = coneSpacingInMeters * 1e6;  % in microns
    
end

function lambda = minConeSpacing(obj)
    mosaicRangeX = obj.center(1) + obj.width/2*[-1 1];
    mosaicRangeY = obj.center(2) + obj.height/2*[-1 1];
    minX = min(abs(mosaicRangeX));
    if (prod(mosaicRangeX) < 0)
        minX = 0;
    end
    minY = min(abs(mosaicRangeY));
    if (prod(mosaicRangeY) < 0)
        minY = 0;
    end
    eccentricityInMeters = sqrt(minX^2 + minY^2);
    ang = atan2(minY, minX)/pi*180;
    [coneSpacingInMeters, aperture, density] = coneSize(eccentricityInMeters,ang);
    lambda = coneSpacingInMeters * 1e6;  % in microns
end

function hexLocs = computeHexGrid(rows, cols, lambda)

    extraCols = round(cols/(sqrt(3)/2)) - cols;
    rectXaxis2 = (1:(cols+extraCols));
    [X2,Y2] = meshgrid(rectXaxis2, 1:rows);
    
    X2 = X2 * sqrt(3)/2;
    for iCol = 1:size(Y2,2)
        Y2(:,iCol) = Y2(:,iCol) - mod(iCol-1,2)*0.5;
    end
    
    % To microns
    densityScalar = sqrt(1.1931428);
    X2 = X2 * lambda * densityScalar;
    Y2 = Y2 * lambda * densityScalar;
    marginInConePositions = 0.1;
    indicesToKeep = (X2>=-marginInConePositions) & ...
                    (X2<= cols+marginInConePositions) &...
                    (Y2>=-marginInConePositions) & ...
                    (Y2<= rows+marginInConePositions);             
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
    
    hexLocs = [xHex(:)-mean(xHex(:))   yHex(:)-mean(yHex(:))];
end

function distances = circularDomainFunction(conePositions, center, radius)
    %  points with positive distance will be excluded
    distances = -(radius-sqrt((conePositions(:,1)-center(1)).^2+(conePositions(:,2)-center(2)).^2));
end

function [coneSpacingInMicrons, eccentricitiesInMicrons] = coneSpacingFunction(conePositions)
    eccentricitiesInMicrons = sqrt(sum(conePositions.^2,2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(conePositions(:,2), conePositions(:,1))/pi*180;
    [coneSpacingInMeters, aperture, density] = coneSize(eccentricitiesInMeters,angles);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6; 
end

function hexLocs = computeHexGridNodesOriginal(obj)

    fprintf('Computing hex mosaic corresponding to rect mosaic with %d rows and %d cols\n', obj.rows, obj.cols);
    extraCols = round(obj.cols/(sqrt(3)/2)) - obj.cols;
    rectXaxis2 = (1:(obj.cols+extraCols));
    [X2,Y2] = meshgrid(rectXaxis2, 1:obj.rows);
    
    X2 = X2 * sqrt(3)/2;
    for iCol = 1:size(Y2,2)
        Y2(:,iCol) = Y2(:,iCol) - mod(iCol-1,2)*0.5;
    end
    
    densityScalar = sqrt(2/sqrt(3.0));
    X2 = X2 * densityScalar;
    Y2 = Y2 * densityScalar;
    marginInConePositions = 0.1;
    indicesToKeep = (X2>=-marginInConePositions) & ...
                    (X2<=obj.cols+marginInConePositions) &...
                    (Y2>=-marginInConePositions) & ...
                    (Y2<=obj.rows+marginInConePositions);             
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
   
    xHex = xHex * obj.patternSampleSize(1);
    yHex = yHex * obj.patternSampleSize(2);
    hexLocs = [xHex(:)-mean(xHex(:))   yHex(:)-mean(yHex(:))];
    
    % move to center
    hexLocs(:,1) = hexLocs(:,1) + obj.center(1);
    hexLocs(:,2) = hexLocs(:,2) + obj.center(2);
   
end


function pattern = rectSampledHexPattern(obj)
    fprintf('\nResampling grid. Please wait ... ');
    % Highres grid
    xRectHiRes = (1:obj.cols) * obj.patternSampleSize(1); xRectHiRes = xRectHiRes - mean(xRectHiRes);
    yRectHiRes = (1:obj.rows) * obj.patternSampleSize(2); yRectHiRes = yRectHiRes - mean(yRectHiRes);
    [xx,yy] = meshgrid(xRectHiRes, yRectHiRes); 
    
    % Match spatial extent
    xRectOriginal = (1:size(obj.patternOriginatingRectGrid,2)) * obj.patternSampleSizeOriginatingRectGrid(1); xRectOriginal = xRectOriginal - mean(xRectOriginal);
    yRectOriginal = (1:size(obj.patternOriginatingRectGrid,1)) * obj.patternSampleSizeOriginatingRectGrid(2); yRectOriginal = yRectOriginal - mean(yRectOriginal);
    xRectOriginal = xRectOriginal / max(xRectOriginal) * max(xRectHiRes);
    yRectOriginal = yRectOriginal / max(yRectOriginal) * max(yRectHiRes);
    [xxx,yyy] = meshgrid(xRectOriginal, yRectOriginal);

    % Generate the high res mosaic pattern
    pattern = zeros(numel(yRectHiRes), numel(xRectHiRes))+1; 
    
    % Determine the closest cone type in the originating  grid
    [~,I] = pdist2([xxx(:) yyy(:)], bsxfun(@minus, obj.coneLocsHexGrid, obj.center), 'euclidean', 'Smallest', 1);
    
    % Determine the closest cone location in the high-res grid
    [~,II] = pdist2([xx(:) yy(:)], bsxfun(@minus, obj.coneLocsHexGrid, obj.center), 'euclidean', 'Smallest', 1);
    
    % That's our cone !
    pattern(ind2sub(size(obj.pattern),II)) = obj.patternOriginatingRectGrid(I);
    fprintf('Done !\n');
    
    % Show patterns (only for debugging purposes)
    debugPlots = false;
    if (debugPlots)
        cmap = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
        figure(10); clf;
        imagesc((obj.patternOriginatingRectGrid-1)/3);
        axis 'equal'; axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'

        figure(11); clf;
        imagesc((pattern-1)/3);
        axis 'equal'; axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'
    end
end