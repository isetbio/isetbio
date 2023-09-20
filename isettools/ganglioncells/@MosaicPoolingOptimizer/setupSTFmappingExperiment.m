function [stimParams, thePresentationDisplay] = setupSTFmappingExperiment(inputConeMosaic, ...
        sceneFOVdegs, retinalImageResolutionDegs, stimulusChromaticity)

    wavelengthSupport = inputConeMosaic.wave;

    switch (stimulusChromaticity)
        case 'achromatic'
            coneContrasts = [1 1 1];
            contrast = 0.75;
        case 'Lcone isolating'
            coneContrasts = [1 0 0];
            contrast = 0.12;
        case 'Mcone isolating'
            coneContrasts = [0 1 0];
            contrast = 0.12;
       case 'Scone isolating'
            coneContrasts = [0 0 1];
            contrast = 0.75;
        otherwise
            error('Unknown stimulus chromaticity: ''%s''.', stimulusChromaticity);
    end

    viewingDistanceMeters = 4;

    deltaOri = 15;
    orientationsTested = 0:deltaOri:(180-deltaOri);
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 12 16 20 24 32 48 64];

    if (max(abs(inputConeMosaic.eccentricityDegs)) > 10)
        spatialFrequenciesTested = [0.06 .125 0.25 0.5 1 2 4 6 8 12 16 20 24 32];
    end

    % Generate a presentation display with a desired resolution
    stimulusPixelsNum = round(max(sceneFOVdegs)/retinalImageResolutionDegs);

    % At least 6 samples / period
    maxSF = 1/(2*3*retinalImageResolutionDegs);
    if (max(spatialFrequenciesTested) > maxSF)
        fprintf('Max SF examined (%2.2f c/deg) is too high for this FOV (%2.2f degs) and pixels num (%d). (SFmax: %2.2f c/deg)\n', ...
            max(spatialFrequenciesTested), max(sceneFOVdegs), stimulusPixelsNum, maxSF);
        idx = find(spatialFrequenciesTested <= maxSF);
        spatialFrequenciesTested = spatialFrequenciesTested(idx);
        if (maxSF > max(spatialFrequenciesTested))
            spatialFrequenciesTested(numel(spatialFrequenciesTested)+1) = maxSF;
        end

        fprintf('Will only measure the STF up to %2.2f c/deg.\n', max(spatialFrequenciesTested));
    end

    stimSizeDegs = max(sceneFOVdegs);
    pixelSizeDegs = retinalImageResolutionDegs;

    thePresentationDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            wavelengthSupport, pixelSizeDegs, ...
            viewingDistanceMeters);

    spatialSupportDegs = rfMappingStimulusGenerator.spatialSupport(...
        stimSizeDegs, pixelSizeDegs);
        
    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', contrast, ...
            'orientationsTested', orientationsTested, ...
            'spatialFrequenciesTested', spatialFrequenciesTested, ...
            'orientationDegs', 0, ...
            'spatialFrequencyCPD', spatialFrequenciesTested(1), ...
            'spatialPhaseIncrementDegs', 30, ...
            'pixelSizeDegs', pixelSizeDegs, ...
            'stimSizeDegs', stimSizeDegs, ...
            'spatialMask', ones(numel(spatialSupportDegs), numel(spatialSupportDegs)), ...
            'wavelengthSupport', displayGet(thePresentationDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(thePresentationDisplay, 'viewing distance') ...
            );

    [~, stimParams.spatialPhasesDegs] = ...
        rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
end

