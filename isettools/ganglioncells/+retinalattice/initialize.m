function [rfPositionsMicrons, radiusMicrons] = initialize(fovDegs, whichEye, params, tStart)
    % Regular hex grid with minimal lambda
    rfPositionsMicrons = generateInitialRFpositions(fovDegs*1.33, params.lambdaMinMicrons);
    
    % Probabilistic sampling according to local RF density
    [rfPositionsMicrons, radiusMicrons] = downSampleInitialRFpositions(rfPositionsMicrons, whichEye, params, tStart);
end

function [rfPositions, radius] = downSampleInitialRFpositions(rfPositions, whichEye, params,  tStart)
    % Set rng seed
    rng(params.rng);
    
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Started with %2.1f million RFs, time lapsed: %f minutes\n', rfsNum/1000000, toc(tStart)/60);
    else
        fprintf('Started with %2.1f thousand RFs, time lapsed: %f minutes\n', rfsNum/1000, toc(tStart)/60);
    end

    fprintf('Removing cones outside the elliptical domain ...');
    radius = max(abs(rfPositions(:)));

    % Remove cones outside the desired region by applying the domain function
    d = params.domainFunction(rfPositions, radius);
    rfPositions = rfPositions(d < params.borderTolerance, :);
    fprintf('... time lapsed: %f minutes.\n', toc(tStart)/60);

    % Compute spacing at all rfPositions
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Computing separations for %2.1f million nodes ...', rfsNum/1000000);
    else
        fprintf('Computing separations for %2.1f thousand nodes ...', rfsNum/1000);
    end
    [rfSpacingMicrons, eccentricitiesMicrons] = params.rfSpacingExactFunction(rfPositions, whichEye);

    fprintf('... time lapsed: %f minutes.',  toc(tStart)/60);

    % sample probabilistically according to coneSpacingFunction
    fprintf('\nProbabilistic sampling ...');
    densityP = RGCmodels.Watson.convert.spacingToDensity(rfSpacingMicrons);
    densityP = densityP / max(densityP(:));
    
    % Remove cones accordingly
    fixedRFPositionsRadiusInCones = 1;

    keptRFIndices = find(...
        (rand(size(rfPositions, 1), 1) < densityP) | ...
        ((eccentricitiesMicrons < fixedRFPositionsRadiusInCones*params.lambdaMinMicrons)) );

    rfPositions = rfPositions(keptRFIndices, :);
    fprintf(' ... done after %f minutes. Kept %d RFs\n', toc(tStart)/60, numel(keptRFIndices));
end
    

function rfPositions = generateInitialRFpositions(fovDegs, lambdaMinMicrons)
    radiusMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(fovDegs/2);
    rows = 2 * radiusMicrons;
    cols = rows;
    rfPositions = computeHexGrid(rows, cols, lambdaMinMicrons);
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