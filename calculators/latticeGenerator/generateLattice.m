function generateLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 1; %30; 
    
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
    
    % Iterative smoothing params
    % 1. Stop if cones move less than this positional tolerance (x gridParams.lambdaMin) in microns
    iterativeParams.dTolerance = 1.0e-4;
    
    % 2. Stop if we exceed this many iterations
    iterativeParams.maxIterations = 3000;
    
    % 3. Trigger Delayun triangularization if rfmovement exceeds this number
    iterativeParams.maxMovementPercentile = 20;
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    iterativeParams.minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    iterativeParams.maxIterationsBeforeRetriangulation = 30;
    
    % 6. Interval to query user whether he/she wants to terminate
    iterativeParams.queryUserIntervalMinutes = 60*12;
    
        
    % STEP 1. Generate initial RF positions in a regular hex lattice with lambda = min separation
    [rfPositions, lambda] = generateInitialRFpositions(mosaicFOVDegs, neuronalType);
    
    % Visualize lattice
    visualizeLattice(rfPositions);
    
    % STEP 2. Downsample uniform grid
    domain.function = @ellipticalDomainFunction;
    domain.ellipseAxes = [1 1.2247];
    domain.maxEcc = max(abs(rfPositions(:)));
    rfPositions = downSampleInitialRFpositions(rfPositions, lambda, domain, neuronalType, whichEye, theRandomSeed);
    visualizeLattice(rfPositions);
    
    % STEP 3. Generate lookup density tables
    [tabulatedDensity, tabulatedEcc] = generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda,  neuronalType, whichEye);
    
    % STEP 4. Iteratively smooth the lattice grid
    [rfPositions] = iterativelySmoothGrid(rfPositions, tabulatedDensity, tabulatedEcc, iterativeParams, lambda, domain, neuronalType, whichEye);
    
end



function distances = ellipticalDomainFunction(rfPositions, maxEcc, ellipseAxes)
    xx = rfPositions(:, 1);
    yy = rfPositions(:, 2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - maxEcc;
end