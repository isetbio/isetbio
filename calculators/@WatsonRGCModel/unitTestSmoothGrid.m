function unitTestSmoothGrid()

    % Generate or view saved mosaic
    generateNewMosaic = true;
    
    
    % Visualize mosaic and progress
    visualizeProgress = generateNewMosaic;

    % Size of mosaic to generate
    mosaicFOVDegs = 20; %30; 
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 48;
    
    % Set a random seed
    theRandomSeed = 1;
    
    % Termination conditions
    % 1. Stop if cones move less than this positional tolerance (x gridParams.lambdaMin) in microns
    dTolerance = 1.0e-4;
    
    % 2. Stop if we exceed this many iterations
    maxIterations = 3000;
    
    % 3. Trigger Delayun triangularization if rfmovement exceeds this number
    maxMovementPercentile = 30;
    
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    maxIterationsBeforeRetriangulation = 30;
    
    % 6. Interval to query user whether he/she wants to terminate
    queryUserIntervalMinutes = 60*12;
    
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    saveFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_maxMovPrctile%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, maxMovementPercentile));

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
    
    gridParams.whichEye = whichEye;
    gridParams.micronsPerDegree = 300;
    gridParams.ellipseAxes = [1 1.2247];
    gridParams.lambdaMin = 2;
    gridParams.borderTolerance = 0.001 * gridParams.lambdaMin;
    gridParams.dTolerance = gridParams.lambdaMin * dTolerance;
    gridParams.rng = theRandomSeed;
    gridParams.maxMovementPercentile = maxMovementPercentile;
    
   
    if (~generateNewMosaic)
        load(saveFileName, 'rfPositionsHistory','iterationsHistory', 'rfPositionsHistory2','iterationsHistory2', 'maxMovements', 'reTriangulationIterations', 'terminationReason');
        fprintf('Termination reason for this mosaic: %s\n', terminationReason)
        hFig = figure(1); clf;
        
        
        renderFinalMosaic = ~true;
        
        if (renderFinalMosaic)
            set(hFig, 'Position', [10 10 950 600], 'Color', [1 1 1]);
            finalRFPositions = squeeze(rfPositionsHistory(end,:,:));
            rfSpacingInMicrons = gridParams.rfSpacingFunctionFull(double(finalRFPositions),  whichEye);
            
            triangleIndices = delaunayn(double(finalRFPositions));
            
            switch (neuronalType)
                case 'cone'
                    rfIDs = generateConeRFIDs(finalRFPositions);
                case 'mRGC'
                    rfIDs = generateRGCRFIDs(finalRFPositions);
                otherwise
                    error('Unknown neuronal type: ''%s''.',  neuronalType);
            end
            
            ax = subplot('Position', [0.04 0.07 0.95 0.925]);
            
            visualizedRect.size = [0.5 0.3];
            plotTriangularizationGrid = ~true;
            
            videoFileName = strrep(saveFileName, '.mat', '_panVideo');
            videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
            videoOBJ.FrameRate = 30;
            videoOBJ.Quality = 100;
            videoOBJ.open();
    
            deltaX = 0.002;
            stepsNum = 15/deltaX;
            for step = 1:1700
                cla(ax);
                visualizedRect.center = [0.49+(step-1)*deltaX 0];
                deltaX = min([deltaX*1.0015 0.01]);
                renderMosaicWithApertures(ax, finalRFPositions, rfIDs, triangleIndices, rfSpacingInMicrons, ...
                    mosaicFOVDegs, gridParams.micronsPerDegree, visualizedRect, neuronalType, plotTriangularizationGrid);
                videoOBJ.writeVideo(getframe(hFig));
            end
            videoOBJ.close();
            
        else
            set(hFig, 'Position', [10 10 1596 1076]);
            generateMosaicProgressVideo(strrep(saveFileName, 'progress', 'video'), hFig , rfPositionsHistory, iterationsHistory, maxMovements, reTriangulationIterations, gridParams.dTolerance, mosaicFOVDegs, gridParams.micronsPerDegree);
        end
        return;
    end
    
    tStart = tic;

    % Generate initial RF positions and downsample according to the density
    rfPositions = generateInitialRFpositions(mosaicFOVDegs*1.07, gridParams.lambdaMin, gridParams.micronsPerDegree);
    [rfPositions, gridParams] = downSampleInitialRFpositions(rfPositions, gridParams, tStart);


    
    
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Iteration: 0, Adusting %2.1f million nodes, time lapsed: %f minutes\n', size(rfPositions,1)/1000000, toc(tStart)/60);
    else
        fprintf('Iteration: 0, Adusting %2.1f thousand nodes, time lapsed: %f minutes\n', size(rfPositions,1)/1000, toc(tStart)/60);
    end
    
    % Precompure table of spacing across eccentricities
    switch (neuronalType)
        case 'cone'
            [tabulatedEccXYMicrons, tabulatedRFSpacingInMicrons] = ...
                computeTableOfConeSpacings(rfPositions, eccentricitySamplesNum, whichEye);
        case 'mRGC'
            [tabulatedRFSpacingInMicrons, tabulatedEccXYMicrons] = ...
                computeTableOfmRGCRFSpacings(rfPositions, eccentricitySamplesNum, whichEye);
        otherwise
            error('Unknown neuronal type: ''%s''.',  neuronalType);
    end
    
    
    % Do it
    [rfPositions, rfPositionsHistory,iterationsHistory, rfPositionsHistory2, iterationsHistory2, maxMovements, reTriangulationIterations, terminationReason] = ...
        smoothGrid(gridParams, rfPositions,  minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation, maxIterations, queryUserIntervalMinutes, ...
        visualizeProgress,  tabulatedEccXYMicrons, tabulatedRFSpacingInMicrons,  mosaicFOVDegs, tStart);        
    
    % Save results
    save(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iterationsHistory', 'rfPositionsHistory2', 'iterationsHistory2', 'maxMovements', 'reTriangulationIterations', ...
        'terminationReason', '-v7.3');
    fprintf('History saved  in %s\n', saveFileName);
end

function rfIDs = generateConeRFIDs(finalRFPositions)

    rfPercentages = [0.62 0.31 0.07];
    ecc = sqrt(sum(finalRFPositions.^2,2));
    rfNums = size(finalRFPositions,1);
    
    LconesNum = floor(rfPercentages(1)*rfNums);
    MconesNum = floor(rfPercentages(2)*rfNums);
    SconesNum = rfNums -LconesNum - MconesNum;
    SconeSkip = round(rfNums/SconesNum);
    MconeSkip = round(rfNums/MconesNum);
    rfIDs(1:rfNums,1) = 1;
    rfIDs(1:SconeSkip:end,1) = 3;
    rfIDs(1:MconeSkip:end,1) = 2;
    rfIDs(find((ecc<30)&(rfIDs==3))) = 2;
end

function rfIDs = generateRGCRFIDs(finalRFPositions)
    rfNums = size(finalRFPositions,1);
    rfIDs(1:rfNums,1) = nan;
end


function [rfPositions, gridParams] = downSampleInitialRFpositions(rfPositions, gridParams,  tStart)
    
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
        fprintf('Computing separations for %2.1f million nodes ...', rfsNum/1000000);
    else
        fprintf('Computing separations for %2.1f thousand nodes ...', rfsNum/1000);
    end
    rfSeparations = feval(gridParams.rfSpacingFunctionFull, rfPositions, gridParams.whichEye);

    fprintf('... time lapsed: %f minutes.',  toc(tStart)/60);

    fprintf('\nProbabilistic sampling ...');
    normalizedRFSeparations = rfSeparations / gridParams.lambdaMin;
    densityP = WatsonRGCModel.densityFromSpacing(normalizedRFSeparations);

    % Remove cones accordingly
    fixedRFPositionsRadiusInCones = 1;
    radii = sqrt(sum(rfPositions.^2,2));

    keptRFIndices = find(...
        (rand(size(rfPositions, 1), 1) < densityP) | ...
        ((radii < fixedRFPositionsRadiusInCones*gridParams.lambdaMin)) );

    rfPositions = rfPositions(keptRFIndices, :);
    fprintf(' ... done ! After %f minutes.\n', toc(tStart)/60);
end
    
function rfPositions = generateInitialRFpositions(fovDegs, lambda, micronsPerDeg)
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

function [tabulatedEccXYMicrons, tabulatedConeSpacingInMicrons] = computeTableOfConeSpacings(rfPositions, eccentricitySamplesNum, whichEye)
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
        % true for a rectangular mosaic. For a hex mosaic, spacing = sqrt(2.0/(sqrt(3)*density)).
        correctionFactor = WatsonRGCModel.spacingFromDensity(1);
        tabulatedConeSpacingInMicrons = correctionFactor*tabulatedConeSpacingInMicrons;
end
    
function [rfPositions, rfPositionsHistory, iterationsHistory, rfPositionsHistory2, iterationsHistory2, maxMovements, reTriangulationIterations, terminationReason] = ...
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
    timePrevious = clock;
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
       % positionalDiffsMetric = prctile(positionalDiffs, 99);
        
        % We need to triangulate again if the positionalDiff is above the set tolerance
%         if ((positionalDiffsMetric > gridParams.positionalDiffToleranceForTriangularization))
%             reTriangulationIsNeeded = true;
%             triangularizationTriggerEvent = 'movement > tolerance';
%         end
        
        % We need to triangulate again if the movement in the current iteration was > the average movement in the last 2 iterations 
        if (numel(maxMovements)>3) && (maxMovements(iteration-1) > 1.02*0.5*(maxMovements(iteration-2)+maxMovements(iteration-3)))
            reTriangulationIsNeeded = true;
             triangularizationTriggerEvent = 'movement stopped decreasing';
        end
        
        % We need to triangulate again if we went for maxIterationsToRetriangulate + some more since last triangularization
        if ((abs(lastTriangularizationAtIteration-iteration-1)) > maxIterationsBeforeRetriangulation+min([10 round(iteration/20)]))
            reTriangulationIsNeeded = true;
            triangularizationTriggerEvent = 'maxIterations passed';
        end
        
        % Do not triangulare if we did one less than minIterationsBeforeRetriangulation before
        if ((abs(lastTriangularizationAtIteration-iteration)) < minIterationsBeforeRetriangulation+min([5 round(iteration/50)]))
            reTriangulationIsNeeded = false;
        end
        
        %
        if (iteration==1)
            reTriangulationIsNeeded = true;
            triangularizationTriggerEvent = '1st iteration';
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
        maxMovement = prctile(movementAmplitudes, gridParams.maxMovementPercentile);
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

            % See if we need to query the user about terminating
            timeLapsedMinutes = etime(clock, timePrevious)/60;
            
            if (timeLapsedMinutes > queryUserIntervalMinutes)
                queryUserWhetherToTerminateSoon = true;
            else
                queryUserWhetherToTerminateSoon = false;
            end
            
            fprintf('\t>Triangularization at iteration: %d/%d (%s) - medianMov: %2.6f, tolerance: %2.3f, time lapsed: %2.1f minutes\n', ...
                iteration, gridParams.maxIterations, triangularizationTriggerEvent, maxMovement, gridParams.dTolerance, toc(tStart)/60);
            
            if (isempty(rfPositionsHistory))
                rfPositionsHistory(1,:,:) = single(oldRFPositions);
                iterationsHistory = iteration;
            else
                rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(oldRFPositions), [1 size(oldRFPositions,1) size(oldRFPositions,2)]));
                iterationsHistory = cat(2, iterationsHistory, iteration);
            end
            
            if (visualizeProgress)
                plotMosaic([], rfPositions, triangleIndices, maxMovements, reTriangulationIterations, histogramDiffWidths, histogramData, checkedBins, gridParams.dTolerance, mosaicFOVDegs, gridParams.micronsPerDegree);
            else
                plotMovementSequence([],maxMovements, gridParams.dTolerance)
                plotMeshQuality([],histogramData, checkedBins, iterationsHistory);
            end
        else
            % Store positions at intermediate iterations (every 5 iterations)
            if (mod(iteration,3)==0)
                if (isempty(rfPositionsHistory))
                    rfPositionsHistory2(1,:,:) = single(rfPositions);
                    iterationsHistory2 = iteration;
                else
                    rfPositionsHistory2 = cat(1, rfPositionsHistory2, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
                    iterationsHistory2 = cat(2, iterationsHistory2, iteration);
                end
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
            timePrevious = clock;
        end
        queryUserWhetherToTerminateSoon = false;
        
        
        if (terminateNowDueToReductionInLatticeQuality)
            % Return the last cone positions
            rfPositions = oldRFPositions;
        end
        
        if (~isempty(userRequestTerminationAtIteration)) && (iteration >= userRequestTerminationAtIteration)
            rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(oldRFPositions), [1 size(oldRFPositions,1) size(oldRFPositions,2)]));
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

function [coneSpacingInMicrons, eccentricitiesInMicrons] = coneSpacingFunctionFull(rfPositions, whichEye)
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(rfPositions(:, 2), rfPositions(:, 1)) / pi * 180;
    coneSpacingInMeters = coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles, 'whichEye', whichEye);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6;
    
    % In ConeSizeReadData, spacing is computed as sqrt(1/density). This is
    % true for a rectangular mosaic. For a hex mosaic, spacing = sqrt(2.0/(sqrt(3)*density)).
    correctionFactor = WatsonRGCModel.spacingFromDensity(1);
    coneSpacingInMicrons = coneSpacingInMicrons*correctionFactor;    
end

function [mRGCSpacingInMicrons, tabulatedEccXYMicrons] = computeTableOfmRGCRFSpacings(rfPositions, eccentricitySamplesNum, whichEye)
    % Compute radial sampling vector of retinal positions that we need to compute spacings for
    eccMicrons = WatsonRGCModel.generateSamplingVectorFromScatteredXYPositions(rfPositions, eccentricitySamplesNum);
    eccMicrons = [-fliplr(eccMicrons) 0 eccMicrons];
    [X,Y] = meshgrid(eccMicrons);
    rfPositions = [X(:) Y(:)];
    
    [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunction(rfPositions, eccentricitySamplesNum, whichEye);
    tabulatedEccXYMicrons = rfPositions;
    mRGCSpacingInMicrons = mRGCSpacingInMicrons';
end

function [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunctionFull(rfPositions, whichEye)
    eccentricitySamplesNum = 128;
    [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunction(rfPositions, eccentricitySamplesNum, whichEye);
end

    
% Function to get RGC spacing based on Watson's mRGC-to-cone ratio and ISETBio's cone density.
function [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunction(rfPositions, eccentricitySamplesNum, whichEye)
    
    switch whichEye
        case 'left'
            theView = 'left eye retina';
        case 'right'
            theView = 'right eye retina';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end

    WatsonRGCModelObj = WatsonRGCModel();
    
    % Compute radial sampling vector of retinal positions that we need to compute density for
    eccMM = 1e-3 *  WatsonRGCModel.generateSamplingVectorFromScatteredXYPositions(rfPositions,eccentricitySamplesNum);
    eccMM = [0 eccMM];
    eccDegs = WatsonRGCModel.rhoMMsToDegs(eccMM);
    
    [coneDensity2DMap, coneMeridianDensities, densitySupportMM, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DConeRFDensity(eccDegs, theView, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    

    [conesToMRGCratio2DMap, spatialSupport, horizontalMeridianLabel, verticalMeridianLabel, ratioLabel, ...
            meridianConeToMRGratio, supportUnits] = WatsonRGCModelObj.compute2DConeToMRGCRFRatio(eccDegs,  theView);
    
    mRGCDensity2DMap = coneDensity2DMap ./ conesToMRGCratio2DMap;
    
    mRGCSpacing2DMapMM = WatsonRGCModelObj.spacingFromDensity(mRGCDensity2DMap);
    
    % Convert to microns from mm
    mRGCSpacing2DMapMicrons = mRGCSpacing2DMapMM*1e3;
    densitySupportMicrons = densitySupportMM*1e3;
    
    % Create a scatterred interpolant function
    [X,Y] = meshgrid(squeeze(densitySupportMicrons(1,:)), squeeze(densitySupportMicrons(2,:)));
    F = scatteredInterpolant(X(:),Y(:),mRGCSpacing2DMapMicrons(:), 'linear');
    
    % Evaluate the interpolant function at the requested rfPositions
    mRGCSpacingInMicrons = F(rfPositions(:,1), rfPositions(:,2));
    
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
end


% Function to get RGC spacing based on Watson's formula. But this is based
% on cone density dropping somewhat differently (in < 0.2 degs) from ISETBio
function [mRGCSpacingInMicrons, eccentricitiesInMicrons] = mRGCSpacingFunctionDirect(rfPositions, eccentricitySamplesNum, whichEye)

    switch whichEye
        case 'left'
            theView = 'left eye retina';
        case 'right'
            theView = 'right eye retina';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end

    WatsonRGCModelObj = WatsonRGCModel();
    
    % Compute radial sampling vector of retinal positions that we need to compute density for
    eccMM = 1e-3 *  WatsonRGCModel.generateSamplingVectorFromScatteredXYPositions(rfPositions,eccentricitySamplesNum);
    eccMM = [0 eccMM];
    
    % Convert retinal mm to visual degs
    eccDegs = WatsonRGCModelObj.rhoMMsToDegs(eccMM);

    % Compute mRGC density map
    [mRGCDensity2DMap, mRGCMeridianDensities, densitySupportMM, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DmRGCRFDensity(eccDegs, theView);


    % Make sure results are returned in retinal mm units
    assert((strcmp(densityUnits, WatsonRGCModelObj.retinalMMDensityUnits)) && ...
           (strcmp(supportUnits, WatsonRGCModelObj.retinalMMEccUnits)), ...
           sprintf('Expected mm units, but got ''%s'' and ''%s'' instead.', supportUnits, densityUnits)); 

    % Density is for both types of mRGCs (ON + OFF), so we need density for
    % one type, which is half (assuming equal numerosities of ON and OFF cells)
    mRGCDensity2DMap = 0.5*mRGCDensity2DMap;
    
    % Convert the density map into a spacing map
    mRGCSpacing2DMapMM = WatsonRGCModelObj.spacingFromDensity(mRGCDensity2DMap);
    
    % Convert to microns from mm
    mRGCSpacing2DMapMicrons = mRGCSpacing2DMapMM*1e3;
    densitySupportMicrons = densitySupportMM*1e3;

    % Create a scatterred interpolant function
    [X,Y] = meshgrid(squeeze(densitySupportMicrons(1,:)), squeeze(densitySupportMicrons(2,:)));
    F = scatteredInterpolant(X(:),Y(:),mRGCSpacing2DMapMicrons(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    mRGCSpacingInMicrons = F(rfPositions(:,1), rfPositions(:,2));
    
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
end


function rfSpacingInMicrons = rfSpacingFunctionFast(rfPositions, tabulatedEccXYMicrons, tabulatedSpacingInMicrons)
    measuresNum = 9;
    [D, I] = pdist2(tabulatedEccXYMicrons, rfPositions, 'euclidean', 'Smallest', measuresNum);

    if (measuresNum > 1)
        totalD = sum(D,1);
        b = zeros(measuresNum, size(D,2));
        for k = 1:measuresNum
            b(k,:) = (totalD - D(k,:)) ./ totalD;
        end
        meanSpacing = sum(b.*tabulatedSpacingInMicrons(I),1); % b1 .* tabulatedSpacingInMicrons(I(1,:)) + b2 .* tabulatedSpacingInMicrons(I(2,:));
    else
        meanSpacing = tabulatedSpacingInMicrons(I);
    end
    
    rfSpacingInMicrons = meanSpacing';
end



function generateMosaicProgressVideo(videoFileName, hFigVideo, rfPositionsHistory, iterationsHistory, maxMovements, reTriangulationIterations, dTolerance, mosaicFOVDegs, micronsPerDegree)
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    widths = [];
    checkedBinHistory = zeros(size(rfPositionsHistory,1), 10);
    for k = 1:size(rfPositionsHistory,1)
        currentRFPositions = squeeze(rfPositionsHistory(k,:,:));
        triangleIndices = delaunayn(double(currentRFPositions));
        [~, histogramData, widths, diffWidths, checkedBins] = checkForEarlyTerminationDueToHexLatticeQualityDecrease(currentRFPositions, triangleIndices, widths);
        checkedBinHistory(k,1:numel(checkedBins)) = checkedBins;
        plotMosaic(hFigVideo, currentRFPositions, triangleIndices, ...
            maxMovements(1:iterationsHistory(k)), reTriangulationIterations(1:k), ...
            diffWidths, histogramData, checkedBins, dTolerance, mosaicFOVDegs, ...
            micronsPerDegree);
        % Add video frame
        videoOBJ.writeVideo(getframe(hFigVideo));
    end
    videoOBJ.close();
    
    figure(12345);
    plot(reTriangulationIterations, checkedBinHistory(:,1), 'rs-');
    xlabel('iteration');
    ylabel('histogram width');
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
    cond2 = (any(diffWidths(:) > 0.1)) && (~any((isnan(diffWidths(:)))));
    if (cond1 && cond2)
        fprintf(2,'Should terminate here\n');
        terminateNow = true;
    else
        terminateNow = false;
    end
        
end

function plotMeshQuality(figNo,histogramData, bin1Percent, iterationsHistory)
    if (isempty(figNo))
        hFig = figure(10); 
        subplotIndex = mod(numel(iterationsHistory)-1,12)+1;
        if (subplotIndex == 1)
            clf;
            set(hFig, 'Position', [10 10 820 930]);
        end
        
        rows = 5; cols = 3;
        posVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rows, ...
           'colsNum', cols, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02);
        row = 1+floor((subplotIndex-1)/cols);
        col = 1+mod((subplotIndex-1),cols);
        subplot('Position', posVectors(row,col).v);
    end
 
    qLims = [0.6 1.005]; 
    bar(histogramData.x,histogramData.y,1); hold on;
    plot(bin1Percent(1)*[1 1], [0 max(histogramData.y)], 'r-', 'LineWidth', 1.5);
    plot(bin1Percent(end)*[1 1], [0 max(histogramData.y)], 'c-', 'LineWidth', 1.5);
    plot(bin1Percent(2)*[1 1], [0 max(histogramData.y)], 'k-',  'LineWidth', 1.5);
    plot(bin1Percent(3)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    plot(bin1Percent(4)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', qLims, 'YLim', [0 max(histogramData.y)], ...
        'XTick', [0.6:0.05:1.0],  'XTickLabel', {'.6', '', '.7', '', '.8', '', '.9', '', '1.'}, ...
        'YTickLabel', {}, 'FontSize', 12);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 12);
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

function renderMosaicWithApertures(axesHandle, rfPositions, rfIDs, triangleIndices, rfSpacing, mosaicFOVDegs, micronsPerDeg, visualizedRect, neuronalType, plotTriangularizationGrid)
    xLimsDeg = visualizedRect.size(1)/2*[-1 1] + visualizedRect.center(1);
    yLimsDeg = visualizedRect.size(2)/2*[-1 1] + visualizedRect.center(2);
    
    
    rfPositionsDeg = rfPositions / micronsPerDeg;
    rfPositionsDegFull = rfPositionsDeg;
    
    rfSpacingDeg = rfSpacing / micronsPerDeg;
    
    idx = find((abs(rfPositionsDeg(:,1)-visualizedRect.center(1))<visualizedRect.size(1)/2) & ...
          (abs(rfPositionsDeg(:,2)-visualizedRect.center(2))<visualizedRect.size(2)/2));
      
    rfPositions = rfPositions(idx,:);
    rfPositionsDeg = rfPositionsDeg(idx,:);
    rfSpacingDeg = rfSpacingDeg(idx);
    rfIDs = rfIDs(idx);
    
    
    % Outline
    switch (neuronalType)
        case 'cone'
             deltaTheta = 20;
             rfOutlineSamplesNum = 360/deltaTheta;
             anglesDegs = (0:rfOutlineSamplesNum)*deltaTheta;
             rDegs = 0.7*(0.5*rfSpacingDeg);
        case 'mRGC'
            deltaTheta = 20;
            rfOutlineSamplesNum = 360/deltaTheta;
            anglesDegs = (0:rfOutlineSamplesNum)*deltaTheta;
            rDegs = 1.0*(0.5*rfSpacingDeg);
    end
    
   
    
    % rf.xProfiles in [nodes x (segmentsNum+1)]
    rf.xProfilesDeg = rDegs * cosd(anglesDegs);
    rf.yProfilesDeg = rDegs * sind(anglesDegs);
   
    
    % Translate profiles according to RF center position
    x = bsxfun(@plus, rf.xProfilesDeg, rfPositionsDeg(:,1));
    y = bsxfun(@plus, rf.yProfilesDeg, rfPositionsDeg(:,2));

    rfContours = cell(1,size(x,1));
    for coneIndex = 1:size(x,1)
       rfContours{coneIndex} = struct('x', x(coneIndex,:), 'y', y(coneIndex,:));  
    end
    
    
    if (plotTriangularizationGrid)
        visualizeLatticeState(rfPositionsDegFull, triangleIndices);
    end
        
    hold on;
    
    edgeColor = [0 0 0]; lineWidth = 1.0;
    renderPatchArray(axesHandle, rfContours, rfIDs, edgeColor, lineWidth);
    
    axis 'equal'
    if (xLimsDeg(1) < 2)
        xTicks = 0:0.1:35;
        yTicks = [-3:0.1:3];
    elseif (xLimsDeg(1) < 4)
        xTicks = 0:0.2:35;
        yTicks = [-3:0.2:3];
    elseif (xLimsDeg(1) <10)
        xTicks = 0:0.5:35;
        yTicks = [-3:0.5:3];
    else
        xTicks = 0:1:35;
        yTicks = [-3:1:3];
    end

    xTicks = 0:0.1:35;
    yTicks = [-3:0.1:3];
        
    box on;
    set(gca, 'Color', 'none', 'XLim', xLimsDeg, 'YLim', yLimsDeg, 'CLim', [1 3], 'FontSize', 16);
    set(gca, 'XTick', xTicks, 'YTick', yTicks);
    xlabel('\it space (degs)');
    drawnow;
    
end

function renderPatchArray(axesHandle, rfContours, rfIDs, edgeColor, lineWidth)
   
    conesNum = numel(rfContours);
    maxVerticesPerRFcontour = 50;
   
    bigN = 0;
    for coneIndex = 1:conesNum
        c = rfContours{coneIndex};
        bigN = bigN + numel(c.x);
    end
    
    verticesList = zeros(bigN, 2);
    facesList = nan(conesNum, maxVerticesPerRFcontour);
    colors = zeros(bigN,1);
    
    idx = 0;
    for coneIndex = 1:conesNum
        c = rfContours{coneIndex};
        N = numel(c.x);
        indices = idx + (1:N);
        verticesList(indices, 1) = c.x;
        verticesList(indices, 2) = c.y;
        facesList(coneIndex,1:N) = indices;
        colors(indices, 1) = rfIDs(coneIndex);
        idx = idx + N;
    end
    
    colorMap = [1 0.2 0.6; 0.3 1 0.6; 0.1 0.2 0.8];
    
    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.FaceAlpha = 0.8;
    S.EdgeColor = [0.3 0.3 0.3];
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
    colormap(axesHandle, colorMap); 
end

function plotMosaic(hFig, rfPositions, triangleIndices, maxMovements,  reTriangulationIterations, widths, histogramData, bin1Percent,  dTolerance, mosaicFOVDegs, micronsPerDeg)

    eccDegs = (sqrt(sum(rfPositions.^2, 2)))/micronsPerDeg;
    idx = find(eccDegs <= min([1 mosaicFOVDegs])/2);
    %idx = 1:size(rfPositions,1);
    
    if (isempty(hFig))
        hFig = figure(1);
        set(hFig, 'Position', [10 10 1596 1076]);
    end
    
    clf;
    subplot(2,3,[1 2 4 5]);
    
    if (1==1)
        plotTriangularizationGrid = true;
        if (plotTriangularizationGrid)
            visualizeLatticeState(rfPositions, triangleIndices);
        end

        plot(rfPositions(idx,1), rfPositions(idx,2), 'r.');
        maxPos = max(max(abs(rfPositions(idx,:))));
        set(gca, 'XLim', maxPos*[-1 1], 'YLim', maxPos*[-1 1], 'FontSize', 16);
        axis 'square'
    end
      
    subplot(2,3,3);
    yyaxis left
    plotMovementSequence(hFig, maxMovements, dTolerance);
    title(sprintf('Iteration: %d', reTriangulationIterations(end)));
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
    
    patch(xx, yy, [0 0 1], 'EdgeColor', [1 0 0], ...
        'EdgeAlpha', 0.8, 'FaceAlpha', 0.0, ...
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
