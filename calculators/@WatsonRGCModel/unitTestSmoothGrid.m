function unitTestSmoothGrid()

    gridParams.coneSpacingFunction = @coneSpacingFunction;
    gridParams.coneSpacingFunctionNew = @coneSpacingFunctionNew;
    gridParams.domainFunction = @ellipticalDomainFunction;
    
    gridParams.center = [0 0];
    gridParams.ellipseAxes = [1 1.2247];
    gridParams.borderTolerance = 0.001 * 2;
    gridParams.positionalDiffTolerance = 0.8;
    gridParams.lambdaMin = 2;
    
    coneLocsDir = '/Users/nicolas/Documents/MATLAB/projects/ISETBioCSF/sideprojects/MosaicGenerator';
    loadConePositions = ~true;
    
    tic
    if (loadConePositions)
        %load('cp0.5degs.mat', 'conePositions');
        
        mosaicFOVDegs = 7.0; % choose from '0.5', '1.5', '7.0';
        load(fullfile(coneLocsDir,sprintf('cp%2.1fdegs.mat',mosaicFOVDegs)), 'conePositions');

        fovMax = 5.0;
        ecc = sqrt(sum(conePositions.^2,2))/300;
        idx = find(ecc <= fovMax/2);
        conePositions = conePositions(idx,:);
        gridParams.radius = max(abs(conePositions(:)));
    else
        mosaicFOVDegs  = 20.0;
        conePositions = generateConePositions(mosaicFOVDegs);
        conesNum = size(conePositions,1);
        if (conesNum > 1000*1000)
            fprintf('Started with %2.1f million cones, time lapsed: %f minutes\n', conesNum/1000000, toc/60);
        else
            fprintf('Started with %2.1f thousand cones, time lapsed: %f minutes\n', conesNum/1000, toc/60);
        end
        
        fprintf('Removing cones outside the ellipse ...');
        gridParams.radius = max(abs(conePositions(:)));
    
        % Remove cones outside the desired region by applying the provided
        % domain function
        d = feval(gridParams.domainFunction, conePositions, ...
            gridParams.center, gridParams.radius, gridParams.ellipseAxes);
        conePositions = conePositions(d < gridParams.borderTolerance, :);
        fprintf('... time lapsed: %f minutes.\n', toc/60);
        
        
        % sample probabilistically according to coneSpacingFunction
        conesNum = size(conePositions,1);
        if (conesNum > 1000*1000)
            fprintf('Computing separations for %2.1f million cones ...', conesNum/1000000);
        else
            fprintf('Computing separations for %2.1f thousand cones ...', conesNum/1000);
        end
        coneSeparations = feval(gridParams.coneSpacingFunction, conePositions);
        fprintf('... time lapsed: %f minutes.',  toc/60);
    
        fprintf('\nProbabilistic sampling ...');
        normalizedConeSeparations = coneSeparations / gridParams.lambdaMin;
        densityP = 1/(sqrt(2/3)) * (1 ./ normalizedConeSeparations) .^ 2;
    
        % Remove cones accordingly
        fixedConePositionsRadiusInCones = 1;
        radii = sqrt(sum(conePositions.^2,2));
    
        keptConeIndices = find(...
            (rand(size(conePositions, 1), 1) < densityP) | ...
            ((radii < fixedConePositionsRadiusInCones*gridParams.lambdaMin)) );
 
        conePositions = conePositions(keptConeIndices, :);
        fprintf(' ... done !\n');
    end
    
    conesNum = size(conePositions,1);
    if (conesNum > 1000*1000)
        fprintf('Iteration: 0, Adusting %2.1f million cones, time lapsed: %f minutes\n', size(conePositions,1)/1000000, toc/60);
    else
        fprintf('Iteration: 0, Adusting %2.1f thousand cones, time lapsed: %f minutes\n', size(conePositions,1)/1000, toc/60);
    end
    
    
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of conePositions
    eccentricitySamplesNum = 32;
    eccSpacePartitions = 4;
    whichEye = 'right';
    
    % Termination conditions
    dTolerance = 1.0e-4;
    maxIterations = 1000;
    
    % Options
    visualizeProgress = ~true;
    useOldMethod = ~true;
    
    % Save filename
    saveFileName = fullfile(coneLocsDir, sprintf('progress_%s_partitionsNum%d_samplesNum_%d.mat', whichEye, eccSpacePartitions, eccentricitySamplesNum));
    

    [tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons] = ...
            computeTableOfConeSpacings(conePositions, eccentricitySamplesNum, whichEye, eccSpacePartitions);
    

    
    % Do it
    [conePositions, conePositionsHistory,iterationsHistory] = ...
        smoothGrid(gridParams, conePositions,  dTolerance, maxIterations, visualizeProgress, eccSpacePartitions, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons, useOldMethod, mosaicFOVDegs);        
    
    % Save results
    save(saveFileName, 'conePositionsHistory','iterationsHistory', ...
        tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons, ...
        '-v7.3');
    fprintf('History saved  in %s\n', saveFileName);
end

function conePositions = generateConePositions(fovDegs)

    micronsPerDeg = 300;
    radius = fovDegs/2*1.2*micronsPerDeg;
    lambda = 2;
    rows = 2 * radius;
    cols = rows;
    conePositions = computeHexGrid(rows, cols, lambda);
end

function hexLocs = computeHexGrid(rows, cols, lambda)
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
end

function [tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons] = computeTableOfConeSpacings(conePositions, eccentricitySamplesNum, whichEye, partitions)
        eccentricitiesInMeters = sqrt(sum(conePositions .^ 2, 2)) * 1e-6;
        s = sort(eccentricitiesInMeters);
        maxConePositionMeters = max(s);
        minConePositionMeters = min(s(s>0));
        eccentricitiesInMeters = logspace(log10(minConePositionMeters), log10(maxConePositionMeters), eccentricitySamplesNum);
        tabulatedEccMeters1D = [-fliplr(eccentricitiesInMeters) 0 eccentricitiesInMeters];
        [tabulatedEccX, tabulatedEccY] = meshgrid(tabulatedEccMeters1D);
        tabulatedEccX = tabulatedEccX(:);
        tabulatedEccY = tabulatedEccY(:);
        
        if (partitions == 1)
            tabulatedEccMeters = sqrt(tabulatedEccX.^2 + tabulatedEccY.^2);
            tabulatedEccAngles = atan2d(tabulatedEccY, tabulatedEccX);
            tabulatedConeSpacingInMeters = coneSizeReadData(...
                'eccentricity', tabulatedEccMeters, ...
                'angle', tabulatedEccAngles, ...
                'whichEye', whichEye);
            
            tabulatedEccXYMicrons = [tabulatedEccX tabulatedEccY]*1e6;
            tabulatedConeSpacingInMicrons = tabulatedConeSpacingInMeters * 1e6;
        else
            for qIndex = 1:4
                switch (qIndex)
                    case 1
                        idx = find( (tabulatedEccX>=0) & (tabulatedEccY>=0));
                    case 2
                        idx = find( (tabulatedEccX>=0) & (tabulatedEccY<=0));
                    case 3
                        idx = find( (tabulatedEccX<=0) & (tabulatedEccY>=0));
                    case 4
                        idx = find( (tabulatedEccX<=0) & (tabulatedEccY<=0));
                end
                
                tabulatedEccMeters = sqrt(tabulatedEccX(idx).^2 + tabulatedEccY(idx).^2);
                tabulatedEccAngles = atan2d(tabulatedEccY(idx), tabulatedEccX(idx));
                tabulatedConeSpacingInMeters = coneSizeReadData(...
                    'eccentricity', tabulatedEccMeters, ...
                    'angle', tabulatedEccAngles, ...
                    'whichEye', whichEye);
                
                tabulatedEccXYMicrons(qIndex,:,:) = [tabulatedEccX(idx) tabulatedEccY(idx)]*1e6;
                tabulatedConeSpacingInMicrons(qIndex,:) = tabulatedConeSpacingInMeters * 1e6;
            end
        end
        
        % In ConeSizeReadData, spacing is computed as sqrt(1/density). This is
        % true for a rectangular mosaic. For a hex mosaic, spacing = sqrt(2.0/(3*density)).
        tabulatedConeSpacingInMicrons = sqrt(2/3)*tabulatedConeSpacingInMicrons;
end
    
function [conePositions, conePositionsHistory,iterationsHistory] = smoothGrid(gridParams, conePositions,  dTolerance, maxIterations, visualizeProgress, eccSpacePartitions, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons, useOldMethod, mosaicFOVDegs)  

    gridParams.maxIterations = maxIterations;
    gridParams.dTolerance = gridParams.lambdaMin * dTolerance;
    deps = sqrt(eps) * gridParams.lambdaMin;
    deltaT = 0.2;

    % Initialize convergence
    forceMagnitudes = [];

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    
    % Number of cones
    conesNum = size(conePositions, 1);
    
    % Iteratively adjust the cone positions until the forces between nodes
    % (conePositions) reach equilibrium.
    notConverged = true;
    oldConePositions = inf;
    
    iteration = 0;
    maxMovements = [];
    conePositionsHistory = [];
    
    tic
    
    while (notConverged) && (iteration <= gridParams.maxIterations)
        iteration = iteration + 1;

        % compute cone positional diffs
        positionalDiffs = sqrt(sum((conePositions-oldConePositions).^ 2,2)); 
        
        reTriangulationIsNeeded = (max(positionalDiffs) > gridParams.positionalDiffTolerance);
        if (reTriangulationIsNeeded)
            % save old come positions
            oldConePositions = conePositions;
            
            % Perform new Delaunay triangulation to determine the updated
            % topology of the truss.
            triangleConeIndices = delaunayn(conePositions);
            % Compute the centroids of all triangles
            centroidPositions = 1.0/3.0 * (...
                    conePositions(triangleConeIndices(:, 1), :) + ...
                    conePositions(triangleConeIndices(:, 2), :) + ...
                    conePositions(triangleConeIndices(:, 3), :));
            
            % Remove centroids outside the desired region by applying the
            % signed distance function
            d = feval(gridParams.domainFunction, centroidPositions, ...
                    gridParams.center, gridParams.radius, ...
                    gridParams.ellipseAxes);
            triangleConeIndices = triangleConeIndices(d < gridParams.borderTolerance, :);
            
           % Create a list of the unique springs (each spring connecting 2 cones)
           springs = [...
                    triangleConeIndices(:, [1, 2]); ...
                    triangleConeIndices(:, [1, 3]); ...
                    triangleConeIndices(:, [2, 3]) ...
           ];
           springs = unique(sort(springs, 2), 'rows');
            
           % find all springs connected to this cone
           springIndices = cell(1,conesNum);
           for coneIndex = 1:conesNum
               springIndices{coneIndex} = find((springs(:, 1) == coneIndex) | (springs(:, 2) == coneIndex));
           end
        end % reTriangulationIsNeeded
        
        % Compute spring vectors
        springVectors =  conePositions(springs(:, 1), :) - conePositions(springs(:, 2), :);
        % their centers
        springCenters = (conePositions(springs(:, 1), :) + conePositions(springs(:, 2), :)) / 2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2, 2));

        
        % Compute desired spring lengths. This is done by evaluating the
        % passed coneDistance function at the spring centers.
        if (useOldMethod)
            desiredSpringLengths = feval(gridParams.coneSpacingFunction, springCenters);
            
        elseif (eccSpacePartitions==1)
            desiredSpringLengths= feval(gridParams.coneSpacingFunctionNew, springCenters, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons);

        elseif (eccSpacePartitions==4)
            desiredSpringLengths = [];
            for qIndex = 1:4
                switch (qIndex)
                    case 1
                        idx = find( (springCenters(:,1)>=0) & (springCenters(:,2)>=0));
                    case 2
                        idx = find( (springCenters(:,1)>=0) & (springCenters(:,2)<=0));
                    case 3
                        idx = find( (springCenters(:,1)<=0) & (springCenters(:,2)>=0));
                    case 4
                        idx = find( (springCenters(:,1)<=0) & (springCenters(:,2)<=0));
                end
                desiredSpringLengths(idx,1) = feval(gridParams.coneSpacingFunctionNew, springCenters(idx,:), squeeze(tabulatedEccXYMicrons(qIndex,:,:)), squeeze(tabulatedConeSpacingInMicrons(qIndex,:)));
            end
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
        netForceVectors = zeros(conesNum, 2);
        
        parfor coneIndex = 1:conesNum
           % compute net force from all connected springs
           deltaPos = -bsxfun(@minus, springCenters(springIndices{coneIndex}, :), conePositions(coneIndex, :));
           netForceVectors(coneIndex, :) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{coneIndex}, :), 1);
        end
            
        % update cone positions according to netForceVectors
        conePositions = conePositions + deltaT * netForceVectors;
        
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
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d < -gridParams.borderTolerance, :) .^2 , 2));
        maxMovement = max(movementAmplitudes);
        maxMovements(iteration) = maxMovement;
        
        if maxMovement < gridParams.dTolerance
            notConverged = false; 
        end
          
        if  ( ...
                (iteration < 10) || ...
                ((iteration < 100)&& (mod(iteration,10) == 0)) || ...
                ((iteration < 500)&& (mod(iteration,50) == 0)) || ...
                ((iteration < 1000)&& (mod(iteration,100) == 0)) || ...
                (mod(iteration,500) == 0) ...
                )
            fprintf('Iteration: %d/%d, maxMov: %2.6f, tolerance: %2.6f, time lapsed: %f minutes\n', ...
                iteration, gridParams.maxIterations, max(movementAmplitudes), gridParams.dTolerance, toc/60);
            
            if (visualizeProgress)
                plotMosaic(conePositions, triangleConeIndices, maxMovements, gridParams.dTolerance, mosaicFOVDegs);
            end
            
            if (isempty(conePositionsHistory))
                conePositionsHistory(1,:,:) = single(conePositions);
                iterationsHistory = iteration;
            else
                conePositionsHistory = cat(1, conePositionsHistory, reshape(single(conePositions), [1 size(conePositions,1) size(conePositions,2)]));
                iterationsHistory = cat(2, iterationsHistory, iteration);
            end
           
        end
    end
    toc
    
    if notConverged
        fprintf('Exceeded max number of iteraritions\n');
    else
        fprintf('Converged !\n');
    end
        
end

function distances = ellipticalDomainFunction(conePositions, center, radius, ellipseAxes)
    xx = conePositions(:, 1) - center(1);
    yy = conePositions(:, 2) - center(2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - radius;
end

function [coneSpacingInMicrons, eccentricitiesInMicrons] = coneSpacingFunction(conePositions)
    eccentricitiesInMicrons = sqrt(sum(conePositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(conePositions(:, 2), conePositions(:, 1)) / pi * 180;
    coneSpacingInMeters = coneSizeReadData('eccentricity', eccentricitiesInMeters, 'angle', angles);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6;
end

function coneSpacingInMicrons = coneSpacingFunctionNew(conePositions, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons)
    [~, I] = pdist2(tabulatedEccXYMicrons, conePositions, 'euclidean', 'Smallest', 1);
    coneSpacingInMicrons = (tabulatedConeSpacingInMicrons(I))';
end

function plotMosaic(conePositions, triangleConeIndices, maxMovements,  dTolerance, mosaicFOVDegs)

    %eccDegs = (sqrt(sum(conePositions.^2, 2)))/300;
    %idx = find(eccDegs <= mosaicFOVDegs/2);
    idx = 1:size(conePositions,1);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1596 1076]);
    subplot(2,3,[1 2 4 5]);
    
    
    plotTriangularizationGrid = true;
    if (plotTriangularizationGrid)
        visualizeLatticeState(conePositions, triangleConeIndices);
    end

    plot(conePositions(idx,1), conePositions(idx,2), 'r.');
    maxPos = max(max(abs(conePositions(idx,:))));
    set(gca, 'XLim', maxPos*[-1 1], 'YLim', maxPos*[-1 1], 'FontSize', 16);
    axis 'square'
    
    subplot(2,3,3);
    if (numel(maxMovements) < 10) 
        markerSize = 12;
    elseif (numel(maxMovements) < 50)
        markerSize = 10;
    elseif (numel(maxMovements) < 100)
        markerSize = 8;
    elseif (numel(maxMovements) < 500)
        markerSize = 6;
    else
        markerSize = 4;
    end
    plot(1:numel(maxMovements), maxMovements, 'ko-', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', markerSize);
    hold on;
    plot([1 numel(maxMovements)], dTolerance*[1 1], 'r-', 'LineWidth', 1.5);
    set(gca, 'YLim', [dTolerance*0.5 max(maxMovements)], 'YScale', 'linear', 'FontSize', 16);
    xlabel('iteration');
    ylabel('max movement', 'FontSize', 16)
    ylabel('movement', 'FontSize', 16);
    axis 'square'
    
    qDist = computeQuality(conePositions, triangleConeIndices);
    qLims = [0 1.005]; qBins = [0.0:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    subplot(2,3,6)
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 max(counts)], 'XTick', [0.1:0.2:1.0],  'FontSize', 16);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('count', 'FontSize', 16);
    axis 'square'
    drawnow
end

function q = computeQuality(coneLocs, triangles)
    
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


function visualizeLatticeState(conePositions, triangleConeIndices)
    x = conePositions(:,1);
    y = conePositions(:,2);
    
    xx = []; yy = [];
    for triangleIndex = 1:size(triangleConeIndices, 1)
        coneIndices = triangleConeIndices(triangleIndex, :);
        xCoords = x(coneIndices);
        yCoords = y(coneIndices);
        for k = 1:numel(coneIndices)
            xx = cat(2, xx, xCoords);
            yy = cat(2, yy, yCoords);
        end
    end
    
    patch(xx, yy, [0 0 1], 'EdgeColor', [0.4 0.4 0.4], ...
        'EdgeAlpha', 0.5, 'FaceAlpha', 0.4, ...
        'FaceColor', [0.99 0.99 0.99], 'LineWidth', 1.0, ...
        'LineStyle', '-', 'Parent', gca); 
    hold on;
end
