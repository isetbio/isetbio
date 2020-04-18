function generateLattice

    % Size of mosaic to generate
    mosaicFOVDegs = 0.5; %30; 
    
    % Type of mosaic to generate
    neuronalType = 'cone';
    neuronalType = 'mRGC';
    
    % Which eye
    whichEye = 'right';
    
    % Samples of eccentricities to tabulate spacing on
    % Precompute cone spacing for a grid of [eccentricitySamplesNum x eccentricitySamplesNum] covering the range of rfPositions
    eccentricitySamplesNum = 32;
    
    % Set a random seed
    theRandomSeed = 1;
    
    % Termination conditions
    % 1. Stop if cones move less than this positional tolerance (x gridParams.lambdaMin) in microns
    dTolerance = 1.0e-4;
    
    % 2. Stop if we exceed this many iterations
    maxIterations = 3000;
    
    % 3. Trigger Delayun triangularization if rfmovement exceeds this number
    maxMovementPercentile = 20;
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    maxIterationsBeforeRetriangulation = 30;
    
    % 6. Interval to query user whether he/she wants to terminate
    queryUserIntervalMinutes = 60*12;
    
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0.1 0.1 0.1; 1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 14, ...
            'figureHeightInches', 14);
        
    % Generate initial RF positions in a regular hex lattice with lambda = min separation
    rfPositions = generateInitialRFpositions(mosaicFOVDegs, neuronalType);
    
    % Visualize lattice
    visualizeLattice(rfPositions);
end

