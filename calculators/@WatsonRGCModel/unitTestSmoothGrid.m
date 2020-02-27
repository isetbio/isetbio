function unitTestSmoothGrid()

    % Generate or view saved mosaic
    generateNewMosaic = true;
    
    % Visualize mosaic and progress
    visualizeProgress = ~generateNewMosaic;

    % Size of mosaic to generate
    mosaicFOVDegs = 15; 
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    %neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 32;
    
    % Set a random seed
    theRandomSeed = 1;
    
    % Termination conditions
    % 1. Stop if cones move less than this positional tolerance (x gridParams.lambdaMin) in microns
    dTolerance = 1.0e-3;
    
    % 2. Stop if we exceed this many iterations
    maxIterations = 3000;
    
    % 3. Trigger Delayun triangularization if rfmovement exceeds this number
    ercentageRFSeparationThresholdForTriangularizationPositionalThreshold = 99;
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    maxIterationsBeforeRetriangulation = 15;
    
    % 6. Interval to query user whether he/she wants to terminate
    queryUserIntervalMinutes = 600;
    
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    saveFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_prctile%d.mat', whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, ercentageRFSeparationThresholdForTriangularizationPositionalThreshold));

    % Set grid params
    switch (neuronalType)
        case 'cone'
            gridParams.rfSpacingFunctionFull = @coneSpacingFunctionFull;
            gridParams.rfSpacingFunctionFast = @rfSpacingFunctionFast;
            gridParams.domainFunction = @ellipticalDomainFunction;
        case 'mRGC'
            gridParams.rfSpacingFunctionFull = @mRGCSpacingFunctionFull;
            gridParams.rfSpacingFunctionFast = @rfSpacingFunctionFast;
            gridParams.domainFunction = @ellipticalDomainFunction;
        otherwise
            error('Unknown neuronal type: ''%s''.',  neuronalType);
    end
    
   
    gridParams.ellipseAxes = [1 1.2247];
    gridParams.lambdaMin = 2;
    gridParams.borderTolerance = 0.001 * gridParams.lambdaMin;
    gridParams.dTolerance = gridParams.lambdaMin * dTolerance;
    gridParams.rng = theRandomSeed;
    
   
    if (~generateNewMosaic)
        load(saveFileName, 'rfPositionsHistory','iterationsHistory', 'maxMovements', 'reTriangulationIterations', 'terminationReason');
        fprintf('Termination reason for this mosaic: %s\n', terminationReason)
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1596 1076]);
        generateMosaicProgressVideo(strrep(saveFileName, 'progress', 'video'), hFig , rfPositionsHistory, iterationsHistory, maxMovements, reTriangulationIterations, gridParams.dTolerance, mosaicFOVDegs);
        return;
    end
    
    % Generate initial RF positions and downsample according to the density
    tStart = tic;
    rfPositions = generateInitialRFpositions(mosaicFOVDegs*1.07, gridParams.lambdaMin);
    [rfPositions, gridParams] = downSampleInitialRFpositions(rfPositions, gridParams, ercentageRFSeparationThresholdForTriangularization, tStart);
       
    
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Iteration: 0, Adusting %2.1f million cones, time lapsed: %f minutes\n', size(rfPositions,1)/1000000, toc(tStart)/60);
    else
        fprintf('Iteration: 0, Adusting %2.1f thousand cones, time lapsed: %f minutes\n', size(rfPositions,1)/1000, toc(tStart)/60);
    end
    
    % Tabulate ecc
    [tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons] = ...
            computeTableOfRFSpacings(rfPositions, eccentricitySamplesNum, whichEye);
    
    % Do it
    [rfPositions, rfPositionsHistory,iterationsHistory, maxMovements, reTriangulationIterations, terminationReason] = ...
        smoothGrid(gridParams, rfPositions,  minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation, maxIterations, queryUserIntervalMinutes, ...
        visualizeProgress,  tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons,  mosaicFOVDegs, tStart);        
    
    % Save results
    save(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iterationsHistory', 'maxMovements', 'reTriangulationIterations', ...
        'terminationReason', 'tabulatedEccXYMicrons', 'tabulatedConeSpacingInMicrons', ...
        '-v7.3');
    fprintf('History saved  in %s\n', saveFileName);
end

function [rfPositions, gridParams] = downSampleInitialRFpositions(rfPositions, gridParams, percentageRFSeparationThresholdForTriangularization, tStart)
    
    rng(gridParams.rng);
    
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Started with %2.1f million RFs, time lapsed: %f minutes\n', rfsNum/1000000, toc(tStart)/60);
    else
        fprintf('Started with %2.1f thousand RFs, time lapsed: %f minutes\n', rfsNum/1000, toc(tStart)/60);
    end

    fprintf('Removing cones outside the ellipse ...');
    gridParams.radius = max(abs(rfPositions(:)));

    % Remove cones outside the desired region by applying the provided
    % domain function
    d = feval(gridParams.domainFunction, rfPositions, ...
        gridParams.radius, gridParams.ellipseAxes);
    rfPositions = rfPositions(d < gridParams.borderTolerance, :);
    fprintf('... time lapsed: %f minutes.\n', toc(tStart)/60);


    % sample probabilistically according to coneSpacingFunction
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Computing separations for %2.1f million cones ...', rfsNum/1000000);
    else
        fprintf('Computing separations for %2.1f thousand cones ...', rfsNum/1000);
    end
    coneSeparations = feval(gridParams.rfSpacingFunctionFull, rfPositions);
    gridParams.positionalDiffToleranceForTriangularization = prctile(coneSeparations,percentageRFSeparationThresholdForTriangularization);

    fprintf('... time lapsed: %f minutes.',  toc(tStart)/60);

    fprintf('\nProbabilistic sampling ...');
    normalizedConeSeparations = coneSeparations / gridParams.lambdaMin;
    densityP = 1/(sqrt(2/3)) * (1 ./ normalizedConeSeparations) .^ 2;

    % Remove cones accordingly
    fixedRFPositionsRadiusInCones = 1;
    radii = sqrt(sum(rfPositions.^2,2));

    keptConeIndices = find(...
        (rand(size(rfPositions, 1), 1) < densityP) | ...
        ((radii < fixedRFPositionsRadiusInCones*gridParams.lambdaMin)) );

    rfPositions = rfPositions(keptConeIndices, :);
    fprintf(' ... done ! After %f minutes.\n', toc(tStart)/60);
end
    
function rfPositions = generateInitialRFpositions(fovDegs, lambda)
    micronsPerDeg = 300;
    radius = fovDegs/2*1.2*micronsPerDeg;
    rows = 2 * radius;
    cols = rows;
    rfPositions = computeHexGrid(rows, cols, lambda);
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
    marginInRFPositions = 0.1;
    indicesToKeep = (X2 >= -marginInRFPositions) & ...
                    (X2 <= cols+marginInRFPositions) &...
                    (Y2 >= -marginInRFPositions) & ...
                    (Y2 <= rows+marginInRFPositions);
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
    hexLocs = [xHex(:) - mean(xHex(:)) yHex(:) - mean(yHex(:))];
end

function [tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons] = computeTableOfRFSpacings(rfPositions, eccentricitySamplesNum, whichEye)
        eccentricitiesInMeters = sqrt(sum(rfPositions .^ 2, 2)) * 1e-6;
        s = sort(eccentricitiesInMeters);
        maxConePositionMeters = max(s);
        minConePositionMeters = min(s(s>0));
        eccentricitiesInMeters = logspace(log10(minConePositionMeters), log10(maxConePositionMeters), eccentricitySamplesNum);
        tabulatedEccMeters1D = [-fliplr(eccentricitiesInMeters) 0 eccentricitiesInMeters];
        [tabulatedEccX, tabulatedEccY] = meshgrid(tabulatedEccMeters1D);
        tabulatedEccX = tabulatedEccX(:);
        tabulatedEccY = tabulatedEccY(:);
        
        tabulatedEccMeters = sqrt(tabulatedEccX.^2 + tabulatedEccY.^2);
        tabulatedEccAngles = atan2d(tabulatedEccY, tabulatedEccX);
        tabulatedConeSpacingInMeters = coneSizeReadData(...
            'eccentricity', tabulatedEccMeters, ...
            'angle', tabulatedEccAngles, ...
            'whichEye', whichEye);

        tabulatedEccXYMicrons = [tabulatedEccX tabulatedEccY]*1e6;
        tabulatedConeSpacingInMicrons = tabulatedConeSpacingInMeters * 1e6;
        
        % In ConeSizeReadData, spacing is computed as sqrt(1/density). This is
        % true for a rectangular mosaic. For a hex mosaic, spacing = sqrt(2.0/(3*density)).
        tabulatedConeSpacingInMicrons = sqrt(2/3)*tabulatedConeSpacingInMicrons;
end
    
function [rfPositions, rfPositionsHistory, iterationsHistory, maxMovements, reTriangulationIterations, terminationReason] = ...
    smoothGrid(gridParams, rfPositions,  minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation, maxIterations, queryUserIntervalMinutes, ...
    visualizeProgress, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons, mosaicFOVDegs, tStart)  

    gridParams.maxIterations = maxIterations;
    deps = sqrt(eps) * gridParams.lambdaMin;
    deltaT = 0.2;

    % Initialize convergence
    forceMagnitudes = [];

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    
    % Number of cones
    rfsNum = size(rfPositions, 1);
    
    % Iteratively adjust the cone positions until the forces between nodes
    % (rfPositions) reach equilibrium.
    notConverged = true;
    terminateNowDueToReductionInLatticeQuality = false;
    oldRFPositions = inf;
    
    iteration = 0;
    maxMovements = [];
    rfPositionsHistory = [];
    
    lastTriangularizationAtIteration = iteration;
    minimalIterationsPerformedAfterLastTriangularization = 0;
    histogramWidths = [];
    reTriangulationIterations = [];
    timeLapsedHoursPrevious = [];
    userRequestTerminationAtIteration = [];
    terminateNow = false;
    
    while (~terminateNow) && (~terminateNowDueToReductionInLatticeQuality) && (notConverged) && (iteration <= gridParams.maxIterations) || ...
            ((lastTriangularizationAtIteration > iteration-minimalIterationsPerformedAfterLastTriangularization)&&(iteration > gridParams.maxIterations))
        
        if ((lastTriangularizationAtIteration > iteration-minimalIterationsPerformedAfterLastTriangularization)&&(iteration > gridParams.maxIterations))
            fprintf('Exceed max iterations (%d), but last triangularization was less than %d iterations before so we will do one more iteration\n', gridParams.maxIterations,minimalIterationsPerformedAfterLastTriangularization);
        end
        
        iteration = iteration + 1;

        % compute cone positional diffs
        positionalDiffs = sqrt(sum((rfPositions-oldRFPositions).^ 2,2)); 
        
        % Check if we need to re-triangulate
        %positionalDiffsMetric = max(positionalDiffs);
        %positionalDiffsMetric = median(positionalDiffs);
        positionalDiffsMetric = prctile(positionalDiffs, 99);
        
        % We need to triangulate again if the positionalDiff is above the set tolerance
        reTriangulationIsNeeded = (positionalDiffsMetric > gridParams.positionalDiffToleranceForTriangularization);
        
        % We need to triangulate again if the movement in the current iteration was > the average movement in the last 2 iterations 
        if (numel(maxMovements)>3) && (maxMovements(iteration-1) > 0.5*(maxMovements(iteration-2)+maxMovements(iteration-3)))
            reTriangulationIsNeeded = true;
        end
        
        % We need to triangulate again if we went for maxIterationsToRetriangulate + some more since last triangularization
        if ((abs(lastTriangularizationAtIteration-iteration-1)) > maxIterationsBeforeRetriangulation+(round(iteration/10)))
            reTriangulationIsNeeded = true;
        end
        
        % Do not triangulare if we did one less than minIterationsBeforeRetriangulation before
        if ((abs(lastTriangularizationAtIteration-iteration-1)) < minIterationsBeforeRetriangulation)
            reTriangulationIsNeeded = false;
        end
        
        %
        if (iteration==1)
            reTriangulationIsNeeded = true;
        end
        
        if (reTriangulationIsNeeded)
            lastTriangularizationAtIteration = iteration;
            % save old come positions
            oldRFPositions = rfPositions;
            
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
            d = feval(gridParams.domainFunction, centroidPositions, ...
                    gridParams.radius, ...
                    gridParams.ellipseAxes);
            triangleIndices = triangleIndices(d < gridParams.borderTolerance, :);
            
           % Create a list of the unique springs (each spring connecting 2 cones)
           springs = [...
                    triangleIndices(:, [1, 2]); ...
                    triangleIndices(:, [1, 3]); ...
                    triangleIndices(:, [2, 3]) ...
           ];
           springs = unique(sort(springs, 2), 'rows');
            
           % find all springs connected to this cone
           springIndices = cell(1,rfsNum);
           for rfIndex = 1:rfsNum
               springIndices{rfIndex} = find((springs(:, 1) == rfIndex) | (springs(:, 2) == rfIndex));
           end
        end % reTriangulationIsNeeded
        
        % Compute spring vectors
        springVectors =  rfPositions(springs(:, 1), :) - rfPositions(springs(:, 2), :);
        % their centers
        springCenters = (rfPositions(springs(:, 1), :) + rfPositions(springs(:, 2), :)) / 2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2, 2));

        if (reTriangulationIsNeeded)
            % Compute desired spring lengths. This is done by evaluating the
            % passed distance function at the spring centers.
            desiredSpringLengths= feval(gridParams.rfSpacingFunctionFast, springCenters, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons);
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
        netForceVectors = zeros(rfsNum, 2);
        
        parfor rfIndex = 1:rfsNum
           % compute net force from all connected springs
           deltaPos = -bsxfun(@minus, springCenters(springIndices{rfIndex}, :), rfPositions(rfIndex, :));
           netForceVectors(rfIndex, :) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{rfIndex}, :), 1);
        end
            
        % update cone positions according to netForceVectors
        rfPositions = rfPositions + deltaT * netForceVectors;
        
        d = feval(gridParams.domainFunction, rfPositions, ...
                gridParams.radius, gridParams.ellipseAxes);
        outsideBoundaryIndices = d > 0;
            
        % And project them back to the domain
        if (~isempty(outsideBoundaryIndices))
                % Compute numerical gradient along x-positions
                dXgradient = (feval(gridParams.domainFunction, ...
                    [rfPositions(outsideBoundaryIndices, 1) + deps, ...
                    rfPositions(outsideBoundaryIndices, 2)], ...
                    gridParams.radius, ...
                    gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                    deps;
                dYgradient = (feval(gridParams.domainFunction, ...
                    [rfPositions(outsideBoundaryIndices, 1), ...
                    rfPositions(outsideBoundaryIndices, 2)+deps], ...
                    gridParams.radius, ...
                    gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                    deps;

                % Project these points back to boundary
                rfPositions(outsideBoundaryIndices, :) = ...
                    rfPositions(outsideBoundaryIndices, :) - ...
                    [d(outsideBoundaryIndices) .* dXgradient, ...
                    d(outsideBoundaryIndices) .* dYgradient];
        end
            
        % Check if all interior nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d < -gridParams.borderTolerance, :) .^2 , 2));
        maxMovement = prctile(movementAmplitudes, 50);
        maxMovements(iteration) = maxMovement;
        
        if maxMovement < gridParams.dTolerance
            notConverged = false; 
        end
          
        % Check for early termination due to decrease in hex lattice quality
        if (reTriangulationIsNeeded)
            reTriangulationIterations = cat(2,reTriangulationIterations, iteration);
            [terminateNowDueToReductionInLatticeQuality, histogramData, histogramWidths, histogramDiffWidths, checkedBins] = ...
                checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions, triangleIndices, histogramWidths);
        end
        
        if  ( reTriangulationIsNeeded || terminateNowDueToReductionInLatticeQuality)  
            % See if another hour passed and asked the used whether to
            % terminate soon
            timeLapsedMinutes = toc(tStart)/60;
            if (isempty(timeLapsedHoursPrevious))
                timeLapsedHoursPrevious = 0;
            end
            
            timeLapsedHours = floor(timeLapsedMinutes/queryUserIntervalMinutes);
            
            if (timeLapsedHours > timeLapsedHoursPrevious)
                queryUserWhetherToTerminateSoon = true;
            else
                queryUserWhetherToTerminateSoon = false;
            end
            timeLapsedHoursPrevious = timeLapsedHours;
            
            fprintf('\t>Iteration: %d/%d, medianMov: %2.6f, tolerance: %2.3f, time lapsed: %f minutes\n', ...
                iteration, gridParams.maxIterations, maxMovement, gridParams.dTolerance, timeLapsedMinutes);
            
            if (isempty(rfPositionsHistory))
                rfPositionsHistory(1,:,:) = single(rfPositions);
                iterationsHistory = iteration;
            else
                rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
                iterationsHistory = cat(2, iterationsHistory, iteration);
            end
            
            if (visualizeProgress)
                plotMosaic([], rfPositions, triangleIndices, maxMovements, reTriangulationIterations, histogramDiffWidths, histogramData, checkedBins, gridParams.dTolerance, mosaicFOVDegs);
            else
                plotMovementSequence([],maxMovements, gridParams.dTolerance)
                plotMeshQuality([],histogramData, checkedBins, iterationsHistory);
            end
        end
        
        if (queryUserWhetherToTerminateSoon)
            fprintf('Another %d minute period has passed. Terminate soon?', queryUserIntervalMinutes);
            userTermination = GetWithDefault(' If so enter # of iteration to terminate on. Otherwise hit enter to continue', 'continue');
            if (~strcmp(userTermination, 'continue'))
                userRequestTerminationAtIteration = str2double(userTermination);
                if (isnan(userRequestTerminationAtIteration))
                    userRequestTerminationAtIteration = [];
                end
            else
                fprintf('OK, will ask again in %d minutes.', queryUserIntervalMinutes);
            end
        end
        queryUserWhetherToTerminateSoon = false;
        
        
        if (terminateNowDueToReductionInLatticeQuality)
            % Return the last cone positions
            rfPositions = rfPositionsLast;
        else
            % Save last rfPositions
            rfPositionsLast = rfPositions;
        end
        
        if (~isempty(userRequestTerminationAtIteration)) && (iteration >= userRequestTerminationAtIteration)
            rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
            iterationsHistory = cat(2, iterationsHistory, iteration);
            reTriangulationIterations = cat(2,reTriangulationIterations, iteration);
            fprintf('Current iteration: %d, user request stop iteration: %d\n', iteration,userRequestTerminationAtIteration)
            terminateNow = true;
        end
        
    end
    
    if (terminateNow)
            terminationReason = sprintf('User requested termination at iteration %d', userRequestTerminationAtIteration);
    else
        if (notConverged)
            if (terminateNowDueToReductionInLatticeQuality)
                terminationReason = 'Decrease in hex lattice quality.';
            else
                terminationReason = 'Exceeded max number of iterations.';
            end
        else
            terminationReason = 'Converged.';
        end
    end
    
    fprintf('Hex lattice adjustment ended. Reason: %s\n', terminationReason);
end

function distances = ellipticalDomainFunction(rfPositions, radius, ellipseAxes)
    xx = rfPositions(:, 1);
    yy = rfPositions(:, 2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - radius;
end

function [coneSpacingInMicrons, eccentricitiesInMicrons] = coneSpacingFunctionFull(rfPositions)
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(rfPositions(:, 2), rfPositions(:, 1)) / pi * 180;
    coneSpacingInMeters = coneSizeReadData('eccentricity', eccentricitiesInMeters, 'angle', angles);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6;
end

function [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunctionFull(rfPositions)

    % Following needs to be replaced with Watson model  calls
%     eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
%     eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
%     angles = atan2(rfPositions(:, 2), rfPositions(:, 1)) / pi * 180;
%     coneSpacingInMeters = coneSizeReadData('eccentricity', eccentricitiesInMeters, 'angle', angles);
%     coneSpacingInMicrons = coneSpacingInMeters' * 1e6;


    mRGCSpacingInMicrons = [];
    eccentricitiesInMicrons = [];
    
end


function rfSpacingInMicrons = rfSpacingFunctionFast(rfPositions, tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons)
    [~, I] = pdist2(tabulatedEccXYMicrons, rfPositions, 'euclidean', 'Smallest', 1);
    rfSpacingInMicrons = (tabulatedConeSpacingInMicrons(I))';
end

function generateMosaicProgressVideo(videoFileName, hFigVideo, rfPositionsHistory, iterationsHistory, maxMovements, reTriangulationIterations, dTolerance, mosaicFOVDegs)
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    widths = [];
    for k = 1:size(rfPositionsHistory,1)
        currentRFPositions = squeeze(rfPositionsHistory(k,:,:));
        triangleIndices = delaunayn(double(currentRFPositions));
        [~, histogramData, widths, diffWidths, checkedBins] = checkForEarlyTerminationDueToHexLatticeQualityDecrease(currentRFPositions, triangleIndices, widths);
        plotMosaic(hFigVideo, currentRFPositions, triangleIndices, maxMovements(1:iterationsHistory(k)), reTriangulationIterations(1:k), diffWidths, histogramData, checkedBins, dTolerance, mosaicFOVDegs);
        % Add video frame
        videoOBJ.writeVideo(getframe(hFigVideo));
    end
    
end

function [terminateNow, histogramData, widths, diffWidths, bin1Percent] = checkForEarlyTerminationDueToHexLatticeQualityDecrease(currentRFPositions, triangleIndices, widths)
    
    qDist = computeQuality(currentRFPositions, triangleIndices);
    qBins = [0.5:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    bin1Percent = prctile(qDist,[0.8 3 7 15 99.8]);
    [~, idx1] = min(abs(centers-bin1Percent(2)));
    [~, idx2] = min(abs(centers-bin1Percent(3)));
    [~, idx3] = min(abs(centers-bin1Percent(4)));
    [~, idxEnd] = min(abs(centers-bin1Percent(end)));
    if (isempty(widths))
        k = 1;
    else
        k = size(widths,1)+1;
    end
    widths(k,:) = centers(idxEnd)-[centers(idx1) centers(idx2) centers(idx3)];
    if (k == 1)
        diffWidths = nan;
    else
        diffWidths = diff(widths,1)./(widths(end,:));
    end

    histogramData.x = centers;
    histogramData.y = counts;

    % Termination condition
    cond1 = bin1Percent(1) > 0.85;
    cond2 = (any(diffWidths(:) > 0.05)) && (~any((isnan(diffWidths(:)))));
    if (cond1 && cond2)
        fprintf(2,'Should terminate here\n');
        terminateNow = true;
    else
        terminateNow = false;
    end
        
end

function plotMeshQuality(figNo,histogramData, bin1Percent, iterationsHistory)
    if (isempty(figNo))
        figure(10); 
        subplotIndex = mod(numel(iterationsHistory)-1,12)+1;
        if (subplotIndex == 1)
            clf;
        end
        subplot(4,3,subplotIndex);
    end
 
    qLims = [0.6 1.005]; 
    bar(histogramData.x,histogramData.y,1); hold on;
    plot(bin1Percent(1)*[1 1], [0 max(histogramData.y)], 'r-', 'LineWidth', 1.5);
    plot(bin1Percent(end)*[1 1], [0 max(histogramData.y)], 'c-', 'LineWidth', 1.5);
    plot(bin1Percent(2)*[1 1], [0 max(histogramData.y)], 'k-',  'LineWidth', 1.5);
    plot(bin1Percent(3)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    plot(bin1Percent(4)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', qLims, 'YLim', [0 max(histogramData.y)], 'XTick', [0.1:0.05:1.0],  'FontSize', 16);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('count', 'FontSize', 16);
    if (isempty(figNo))
        title(sprintf('iteration:%d', iterationsHistory(end)))
        drawnow;
    end
    
    if (isempty(figNo))
        figure(11); hold on;
    end
    
end

function plotMovementSequence(figNo, maxMovements, dTolerance)
    if (isempty(figNo))
        figure(11); clf;
    end
    
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
    set(gca, 'YLim', [dTolerance*0.5 max(maxMovements)], 'YScale', 'log', 'FontSize', 16);
    xlabel('iteration');
    ylabel('median movement', 'FontSize', 16)
end


function plotMosaic(hFig, rfPositions, triangleIndices, maxMovements,  reTriangulationIterations, widths, histogramData, bin1Percent,  dTolerance, mosaicFOVDegs)

    eccDegs = (sqrt(sum(rfPositions.^2, 2)))/300;
    idx = find(eccDegs <= min([1 mosaicFOVDegs])/2);
    %idx = 1:size(rfPositions,1);
    
    if (isempty(hFig))
        hFig = figure(1);
        set(hFig, 'Position', [10 10 1596 1076]);
    end
    
    clf;
    subplot(2,3,[1 2 4 5]);
    plotTriangularizationGrid = true;
    if (plotTriangularizationGrid)
        visualizeLatticeState(rfPositions, triangleIndices);
    end

    plot(rfPositions(idx,1), rfPositions(idx,2), 'r.');
    maxPos = max(max(abs(rfPositions(idx,:))));
    set(gca, 'XLim', maxPos*[-1 1], 'YLim', maxPos*[-1 1], 'FontSize', 16);
    axis 'square'
    
   
    
    subplot(2,3,3);
    yyaxis left
    plotMovementSequence(hFig, maxMovements, dTolerance);
    
    yyaxis right
    if (~isnan(widths))
        plot(reTriangulationIterations(2:end), widths(:,1), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10); hold on
        plot(reTriangulationIterations(2:end), widths(:,2), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10);
        plot(reTriangulationIterations(2:end), widths(:,3), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10);
    end
    set(gca, 'YLim', [-1.5 0.1]);
     
    subplot(2,3,6);
    plotMeshQuality(hFig,histogramData, bin1Percent, []);
    drawnow
end

function q = computeQuality(rfLocs, triangles)
    
    trianglesNum = size(triangles,1);
    X = rfLocs(:,1);
    Y = rfLocs(:,2);
    
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


function visualizeLatticeState(rfPositions, triangleIndices)
    x = rfPositions(:,1);
    y = rfPositions(:,2);
    
    xx = []; yy = [];
    for triangleIndex = 1:size(triangleIndices, 1)
        coneIndices = triangleIndices(triangleIndex, :);
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
