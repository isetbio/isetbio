function [theScene, theNullScene] = generateSinusoidalStimulusAndBackgroundScene(...
    presentationDisplay, stimParams)


    stimParams.sigmaDegs = stimParams.sizeDegs * 0.15;
    if (~stimParams.isWindowed)
        stimParams.sigmaDegs = 1e9;
    end

    % Remove fields that are not known to generateGaborScene
    stimParams = rmfield(stimParams, 'isWindowed');
    stimParams = rmfield(stimParams, 'positionDegs');

    theScene = generateGaborScene(...
        'stimParams', stimParams, ...
        'presentationDisplay', presentationDisplay);
    
    stimParams.contrast = 0.0;
    theNullScene = generateGaborScene(...
        'stimParams', stimParams, ...
        'presentationDisplay', presentationDisplay);
end