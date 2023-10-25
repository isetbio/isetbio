function [theScene, theNullScene] = generateSinusoidalStimulusAndBackgroundScene(...
    presentationDisplay, sizeDegs, spatialFrequencyCyclesPerDeg, contrast, isWindowed)

    sigmaDegs = sizeDegs * 0.15;
    if (~isWindowed)
        sigmaDegs = 1e6;
    end

    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', spatialFrequencyCyclesPerDeg, ... 
        'orientationDegs', 0, ...               % 45 degrees
        'phaseDegs', 0, ...                     % spatial phase in degrees
        'sizeDegs', sizeDegs, ...               % size in degrees
        'sigmaDegs', sigmaDegs, ...             % sigma of Gaussian envelope
        'contrast', contrast,...                     % 0.6 Michelson contrast
        'meanLuminanceCdPerM2', 40, ...         % mean luminance
        'pixelsAlongWidthDim', 1024, ...         % pixels- width dimension
        'pixelsAlongHeightDim', 1024 ...         % pixels- height dimension
    );

    theScene = generateGaborScene(...
        'stimParams', stimParams, ...
        'presentationDisplay', presentationDisplay);
    
    stimParams.contrast = 0.0;
    theNullScene = generateGaborScene(...
        'stimParams', stimParams, ...
        'presentationDisplay', presentationDisplay);
end