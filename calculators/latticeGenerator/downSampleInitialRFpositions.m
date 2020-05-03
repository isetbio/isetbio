function rfPositions = downSampleInitialRFpositions(rfPositions, lambda, domain, neuronalType, whichEye, theRandomSeed)
    % Remove rfs outside the desired region by applying the provided domain function
    d = feval(domain.function, rfPositions, domain.maxEcc, domain.ellipseAxes);
    rfPositions = rfPositions(d < lambda/1000,:);

    % Determine smallest spacing (delta)
    switch (neuronalType)
        case 'cone'
            [densityPerMM2, maxDensity] = coneDensityFunctionFull(rfPositions, whichEye);
        case 'mRGC'
            [densityPerMM2, maxDensity] = mRGCRFDensityFunctionFull(rfPositions, whichEye);
        otherwise
            error('Unknown neuronalType: ''%s''.', neuronalType)
    end
    
    % Compute normalized densities
    normDensities = densityPerMM2/maxDensity;
    rfsNum = size(rfPositions, 1);
    ecc = sqrt(sum(rfPositions.^2,2));
    
    % Reset the rng
    rng(theRandomSeed);

    % Probabilistically remove rfs with a P(remove) = 1 - norm density
    keptRFIndices = find((rand(rfsNum, 1) < normDensities') | (ecc < lambda)) ;
    rfPositions = rfPositions(keptRFIndices, :);
    fprintf('Kept %2.2f%% of the RFs (n: %d, mean density:%2.3f) \n', size(rfPositions,1)/rfsNum*100, size(rfPositions,1), mean(normDensities));
end

