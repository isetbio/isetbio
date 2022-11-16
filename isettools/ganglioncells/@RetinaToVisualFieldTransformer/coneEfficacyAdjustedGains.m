function [efficacyAdjustedConeWeights, maxEfficacy] = coneEfficacyAdjustedGains(theConeMosaic, ...
           coneDiametersDegs, outerSegmentLengthEccVariationAttenuationFactors, coneWeights, maxEfficacy)

    % Logic from cMosaic.computeAbsorptionRate(), lines 42-55
    if (isfield(theConeMosaic.coneApertureModifiers, 'shape'))
        switch (theConeMosaic.coneApertureModifiers.shape)
            case 'Gaussian'
                characteristicRadiiDegs =  theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor  * coneDiametersDegs;
                apertureAreasDegsSquared = pi * (characteristicRadiiDegs).^2;
            case 'Pillbox'
                apertureAreasDegsSquared = pi * (0.5*coneDiametersDegs).^2;
            otherwise
                error('Do not know how to generate a ''%s'' aperture.', theConeMosaic.coneApertureModifiers.apertureShape)
        end
    else
        % By default we are using the pillbox aperture area
        apertureAreasDegsSquared = pi * (0.5*coneDiametersDegs).^2;
    end

    absorptionEfficacies = apertureAreasDegsSquared .* outerSegmentLengthEccVariationAttenuationFactors;
    if (isempty(maxEfficacy))
        maxEfficacy = max(absorptionEfficacies(:));
    end
    relativeAbsorptionEfficacies = absorptionEfficacies / maxEfficacy;

    % Adjust cone weight by dividing with the relative absorption  efficacy of that cone
    efficacyAdjustedConeWeights = coneWeights(:) ./ relativeAbsorptionEfficacies(:); 
    efficacyAdjustedConeWeights = reshape(efficacyAdjustedConeWeights, size(coneWeights));
end