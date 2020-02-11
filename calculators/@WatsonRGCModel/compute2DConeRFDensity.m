function [coneRFDensity2D, spatialSupport] = compute2DConeRFDensity(obj, eccDegsInREVisualSpace, eccUnits, densityUnits)

    eccSamplesNum = numel(eccDegsInREVisualSpace);
    [eccXInREVisualSpace, eccYInREVisualSpace] = meshgrid(eccDegsInREVisualSpace,eccDegsInREVisualSpace);
    eccXInREVisualSpace = reshape(eccXInREVisualSpace, [eccSamplesNum*eccSamplesNum 1]);
    eccYInREVisualSpace = reshape(eccYInREVisualSpace, [eccSamplesNum*eccSamplesNum 1]);
    
    requestedEccentricities = sqrt(eccXInREVisualSpace.^2 + eccYInREVisualSpace.^2);
    requestedAngles = atan2d(eccYInREVisualSpace,eccXInREVisualSpace);
    
    meridianConeRFSpacing = zeros(numel(obj.enumeratedMeridianNames), numel(requestedEccentricities));
    meridianConeRFDensity = zeros(numel(obj.enumeratedMeridianNames), numel(requestedEccentricities));
    
    % Compute variation along each of the enumerated meridians
    for meridianIndex = 1:numel(obj.enumeratedMeridianNames)
        [meridianConeRFSpacing(meridianIndex,:), meridianConeRFDensity(meridianIndex,:)] = ...
            obj.coneRFSpacingAndDensity(requestedEccentricities, obj.enumeratedMeridianNames{meridianIndex}, eccUnits, densityUnits); 
    end
    
    
    % Do angular interpolation
    coneRFDensity2D = obj.interpolatedValuesFromMeridianValues(meridianConeRFDensity, requestedAngles);
    coneRFDensity2D = reshape(coneRFDensity2D, [eccSamplesNum eccSamplesNum]);
    spatialSupport(1,:) = eccDegsInREVisualSpace;
    spatialSupport(2,:) = eccDegsInREVisualSpace;
end
