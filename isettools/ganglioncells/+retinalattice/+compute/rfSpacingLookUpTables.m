function [tabulatedRFspacingMicrons, tabulatedEccXYMicrons] = ...
    rfSpacingLookUpTables(rfPositionsMicrons, whichEye, rfSpacingExactFunction, lookUpTableSamplesNum)
    
    % Compute radial sampling vector of retinal positions that we need to compute spacings for
    eccMicrons = logSamplingVectorFromScatteredXYPositions(rfPositionsMicrons, lookUpTableSamplesNum);
    eccMicrons = [-fliplr(eccMicrons) 0 eccMicrons];
    [X,Y] = meshgrid(eccMicrons);
    tabulatedEccXYMicrons = [X(:) Y(:)];
    
    rfSpacingMicrons  = rfSpacingExactFunction(tabulatedEccXYMicrons, whichEye);
    tabulatedRFspacingMicrons = rfSpacingMicrons;
end

function samplingVector = logSamplingVectorFromScatteredXYPositions(scatteredXYPositionsMicrons, nSamples)

    eccentricitiesInMicrons = sqrt(sum(scatteredXYPositionsMicrons .^ 2, 2));
    minSeparationMicrons = 2;
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Samping vector
    samplingVector = logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples);
end