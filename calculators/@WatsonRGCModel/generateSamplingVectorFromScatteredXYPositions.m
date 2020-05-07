function samplingVector = generateSamplingVectorFromScatteredXYPositions(scatteredXYPositionsMicrons, nSamples)

    eccentricitiesInMicrons = sqrt(sum(scatteredXYPositionsMicrons .^ 2, 2));
    minSeparationMicrons = 2;
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Samping vector
    samplingVector = logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples);
end


