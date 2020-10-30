function p = configure(fovDegs, neuronType, whichEye)
    
    % Validate input
    validNeuronTypes = retinalattice.validvalues.neuronTypes;
    validEyes = retinalattice.validvalues.eyes; 
    assert(ismember(neuronType, validNeuronTypes), sprintf('Unknown neuron type: ''%s''.', neuronType));
    assert(ismember(whichEye, validEyes), sprintf('Unknown eye: ''%s''.', whichEye));
    
    % Lattice gallery directory
    p.latticeGalleryDir = sprintf('%s/%s/%s', isetbioRootPath, 'isettools/ganglioncells/latticegallery');

    % Patch filename
    p.patchSaveFileName = sprintf('%s_%s_%1.0fdeg_mosaic_progress', ...
        strrep(whichEye, ' ', '_'), strrep(neuronType, ' ', '_'), fovDegs);
    
    % Radnom number generator seed
    p.rng = 1;
    
    % Max number of iterations
    p.maxIterations = 3*1000;
    
    % Save positions every 5 iterations
    p.iterationIntervalForSavingPositions = 5;
    
    % Terminate iterative method if we achieve at least this mesh quality level
    p.minHexQualityForTermination = 0.84;
    
    % Percetile of q-distribution that we track for convergence
    p.qDistPercentile = 0.8;
    
    % How often to query user whether to terminate iterations
    % p.queryUserIntervalMinutes = 60*12;
    
    % Tolerance for moving point back to elliptical domain
    p.lambdaMinMicrons = 2;
    p.borderTolerance = 0.001 * p.lambdaMinMicrons;
    
    % Terminate loop if rfs move less than this positional tolerance 
    p.dTolerance = 1.0e-4 * p.lambdaMinMicrons;
    
    % Do not trigger Delayun triangularization if < 5 iterations since last one
    p.minIterationsBeforeRetriangulation = 5;
    
    % Trigger Delayun triangularization if > 30 iterations since last one
    p.maxIterationsBeforeRetriangulation = 15;
    
    % Trigger Delayun triangularization if rfmovement exceeds this tolerance
    p.maxMovementPercentile = 30;
    
    % Function handle for domain
    p.domainFunction = @ellipticalDomainFunction;
    
    % Function handle for fast lookup of rf spacing at queried positions
    p.rfSpacingFastFunction = @rfSpacingFromTable;
    
    % Samples for the logarithmically sampled eccentricity lookup table
    p.eccentricityLookUpTableSamplesNum = 48;
    
    % Function handle for computing exact rf spacing at queried positions
    switch (neuronType)
        case 'cones'
            p.rfSpacingExactFunction = @coneSpacingFunction;
        case 'midget ganglion cells'
            p.rfSpacingExactFunction = @midgetRGCSpacingFunction;
    end
end


function [rfSpacingMicrons, eccentricitiesMicrons] = midgetRGCSpacingFunction(rfPosMicrons, whichEye)
    [~, rfSpacingMMs] = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(...
        whichEye, rfPosMicrons, 'midget ganglion cells');
    
    rfSpacingMicrons = 1e3 * rfSpacingMMs(:); 
    eccentricitiesMicrons = sqrt(sum(rfPosMicrons .^ 2, 2));
end

function [rfSpacingMicrons, eccentricitiesMicrons] = coneSpacingFunction(rfPosMicrons, whichEye)
    [~, rfSpacingMMs] = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(...
        whichEye, rfPosMicrons, 'cones');
    
    rfSpacingMicrons = 1e3 * rfSpacingMMs(:); 
    eccentricitiesMicrons = sqrt(sum(rfPosMicrons .^ 2, 2));
end

function distances = ellipticalDomainFunction(rfPositions, radius)
    ellipseAxes = [1 1.2247];
    xx = rfPositions(:, 1);
    yy = rfPositions(:, 2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - radius;
end

function rfSpacingMicrons = rfSpacingFromTable(rfPositions, tabulatedEccXYMicrons, tabulatedSpacingInMicrons)
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
    
    rfSpacingMicrons = meanSpacing';
end
