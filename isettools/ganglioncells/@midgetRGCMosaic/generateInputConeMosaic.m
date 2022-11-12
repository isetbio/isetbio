function generateInputConeMosaic(obj, pResults)
    % No input cone mosaic was passed. We will generate one.
    % Parse inputs relevant to the cone mosaic.
    if (~isempty(pResults.eccentricityDegs))
        % eccentricityDegs is always used if present, for
        % historical reasons
        obj.eccentricityDegs = pResults.eccentricityDegs;
    elseif ~isempty(pResults.positionDegs)
        % Use positionDegs if available but there is not
        % eccentricityDegs
        obj.eccentricityDegs = p.Results.positionDegs;
    else
        % Default
        obj.eccentricityDegs = [0,0];
    end
    
    obj.sizeDegs = pResults.sizeDegs;
    obj.whichEye = pResults.whichEye;
    
    % Set cone aperture modifiers
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
        'smoothLocalVariations', true, ...
        'sigma',  sigmaGaussian, ...
        'shape', 'Gaussian');
    
    % Make the inputConeMosaic a little bit larger than the
    % midget RGC mosaic to allow for cone inputs to the RF surrounds
    obj.extraDegsForInputConeMosaic = midgetRGCMosaic.extraConeMosaicDegsForMidgetRGCSurrounds(...
        obj.eccentricityDegs, obj.sizeDegs);
    
    % Generate the input cone mosaic
    obj.inputConeMosaic = cMosaic(...
        'sourceLatticeSizeDegs', obj.sourceLatticeSizeDegs, ...
        'eccentricityDegs', obj.eccentricityDegs, ...
        'sizeDegs', obj.sizeDegs + obj.extraDegsForInputConeMosaic, ...
        'whichEye', obj.whichEye, ...
        'customDegsToMMsConversionFunction', pResults.customDegsToMMsConversionFunction, ...
        'customMMsToDegsConversionFunction', pResults.customMMsToDegsConversionFunction, ...
        'overlappingConeFractionForElimination', 0.5, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', coneApertureModifiers);
end
