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
    
    % Input cone mosaic generation
    if (isempty(pResults.withInputConeMosaic))
        % Make the inputConeMosaic a little bit larger than the
        % midget RGC mosaic to allow for cone inputs to the RF surrounds
        obj.extraDegsForInputConeMosaic = mRGCMosaic.extraConeMosaicDegsForMidgetRGCSurrounds(...
            obj.eccentricityDegs, obj.sizeDegs);

        % Set cone aperture modifiers
        % Use a Gaussian cone aperture with
        % sigma equal to 0.204 x inner segment diameter (cone diameter)
        sigmaGaussian = 0.204;  % From McMahon et al, 2000
        coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');
        
        % Generate the input cone mosaic
        obj.inputConeMosaic = cMosaic(...
            'sourceLatticeSizeDegs', obj.sourceLatticeSizeDegs, ...
            'eccentricityDegs', obj.eccentricityDegs, ...
            'sizeDegs', obj.sizeDegs + obj.extraDegsForInputConeMosaic, ...
            'whichEye', obj.whichEye, ...
            'eccVaryingConeBlur', true, ...
            'customDegsToMMsConversionFunction', pResults.customDegsToMMsConversionFunction, ...
            'customMMsToDegsConversionFunction', pResults.customMMsToDegsConversionFunction, ...
            'overlappingConeFractionForElimination', 0.5, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', coneApertureModifiers, ...
            'randomSeed', pResults.randomSeed);
    else
        % Use passed input cone mosaic
        obj.inputConeMosaic = pResults.withInputConeMosaic;
        obj.extraDegsForInputConeMosaic = obj.inputConeMosaic.sizeDegs(1) - obj.sizeDegs(1);
        fprintf(2,'Using passed input cone mosaic with extraDegs = %2.2f degs\n', obj.extraDegsForInputConeMosaic);
    end

end