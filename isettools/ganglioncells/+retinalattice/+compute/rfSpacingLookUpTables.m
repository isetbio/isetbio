function [tabulatedRFspacingMicrons, tabulatedEccXYMicrons] = ...
    rfSpacingLookUpTables(rfPositionsMicrons, whichEye, useParfor, rfSpacingExactFunction, lookUpTableSamplesNum)
    
    % Compute radial sampling vector of retinal positions that we need to compute spacings for
    eccMicrons = logSamplingVectorFromScatteredXYPositions(rfPositionsMicrons, lookUpTableSamplesNum);
    %eccMicrons = linSamplingVectorFromScatteredXYPositions(rfPositionsMicrons, lookUpTableSamplesNum);
    eccMicrons = [-fliplr(eccMicrons) 0 eccMicrons];
    [X,Y] = meshgrid(eccMicrons);
    tabulatedEccXYMicrons = [X(:) Y(:)];
    
    rfSpacingMicrons  = rfSpacingExactFunction(tabulatedEccXYMicrons, whichEye, useParfor);
    tabulatedRFspacingMicrons = rfSpacingMicrons;
end

function samplingVector = linSamplingVectorFromScatteredXYPositions(scatteredXYPositionsMicrons, nSamples)

    eccentricitiesInMicrons = sqrt(sum(scatteredXYPositionsMicrons .^ 2, 2));
    minSeparationMicrons = 2;
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Samping vector
    samplingVector = linspace(minSeparationMicrons, maxEccMicrons+minSeparationMicrons, nSamples);
end

function samplingVector = logSamplingVectorFromScatteredXYPositions(scatteredXYPositionsMicrons, nSamples)

    eccentricitiesInMicrons = sqrt(sum(scatteredXYPositionsMicrons .^ 2, 2));
    minSeparationMicrons = 2;
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Samping vector
    samplingVector = logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples);
end