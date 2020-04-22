function generateLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 1.0; %30; 
    visualizedFOVMicrons = 200;
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    %neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 32;
    
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    saveFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum));

    
    % Set a random seed
    theRandomSeed = 1;
    
    % Iterative smoothing params
    % 1. Stop if rfs move less than this positional tolerance in microns
    iterativeParams.dTolerance = 2.0e-4;
    % 1a. Trigger Delayun triangularization if rfmovement exceeds this number
    iterativeParams.maxMovementPercentile = 20;
    
    % 2. Stop if we exceed this many iterations
    iterativeParams.maxIterations = 3000;
    
    %3. Trigger Delayun triangularization if
    %(rfspacing-desiredSpacing)/desiredSpacing > threshold
    iterativeParams.thresholdSpacingDeviation = 0.4;
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    iterativeParams.minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    iterativeParams.maxIterationsBeforeRetriangulation = 30;
    
    % 6. Interval to query user whether he/she wants to terminate
    iterativeParams.queryUserIntervalMinutes = 60*12;
    
    % STEP 1. Generate initial RF positions in a regular hex lattice with lambda = min separation
    [rfPositions, lambda] = generateInitialRFpositions(mosaicFOVDegs, neuronalType);
    
    % STEP 2. Downsample uniform grid
    domain.function = @ellipticalDomainFunction;
    domain.ellipseAxes = [1 1.2247];
    domain.maxEcc = max(abs(rfPositions(:)));
    rfPositions = downSampleInitialRFpositions(rfPositions, lambda, domain, neuronalType, whichEye, theRandomSeed);
    
    % STEP 3. Generate lookup density tables
    [tabulatedDensity, tabulatedSpacing, tabulatedEcc] = generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda,  neuronalType, whichEye);
    
    % STEP 4. Iteratively smooth the lattice grid
    [rfPositions, rfPositionsHistory, maxMovements, iteration] = iterativelySmoothLattice(rfPositions, tabulatedSpacing, tabulatedEcc, iterativeParams, lambda, domain, visualizedFOVMicrons);
    
    % STEP 5. Save results
    save(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iteration', 'maxMovements', '-v7.3');
    fprintf('History saved  in %s\n', saveFileName);
end



function distances = ellipticalDomainFunction(rfPositions, maxEcc, ellipseAxes)
    xx = rfPositions(:, 1);
    yy = rfPositions(:, 2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - maxEcc;
end