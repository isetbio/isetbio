function [coneRFDensity2D, meridianDensities, spatialSupport, xLabelString, yLabelString, ...
    densityLabelString, eccUnits, densityUnits] = compute2DConeRFDensity(obj, eccDegsInREVisualSpace, theReturnedView, varargin)
    
    % Parse input
    p = inputParser;       
    p.addParameter('correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', true, @islogical)
    p.parse(varargin{:});
    correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio = ...
        p.Results.correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio;
    
    % Validate the view name
    obj.validateViewName(theReturnedView);
    
    switch (theReturnedView)
        case {obj.rightEyeVisualField}
            eccUnits = obj.visualDegsEccUnits;
            densityUnits = obj.visualDegsDensityUnits;
            ecc = eccDegsInREVisualSpace;
        case {obj.rightEyeRetina, obj.leftEyeRetina}
            eccUnits = obj.retinalMMEccUnits;
            densityUnits = obj.retinalMMDensityUnits;
            ecc = obj.rhoDegsToMMs(eccDegsInREVisualSpace);
    end
    
    % Symmetrical around 0 eccentricity
    ecc = ecc(2:end);
    ecc = [-fliplr(ecc) 0 ecc];
     
    % Returned density label
    densityLabelString = sprintf('density (cones/%s)',densityUnits);
    
    % Returned spatial support
    spatialSupport(1,:) = ecc;
    spatialSupport(2,:) = ecc;
    
    % Generate (XY) grid of eccentricities
    eccSamplesNum = numel(ecc);
    [eccX, eccY] = meshgrid(ecc,ecc);
    eccX = reshape(eccX, [eccSamplesNum*eccSamplesNum 1]);
    eccY = reshape(eccY, [eccSamplesNum*eccSamplesNum 1]);
    
    % Ecc radii and angles
    requestedEccentricities = sqrt(eccX.^2 + eccY.^2);
    requestedAngles = atan2d(eccY,eccX);
    
    % Preallocate memory
    meridianConeRFSpacing = zeros(numel(obj.enumeratedMeridianNames), numel(requestedEccentricities));
    meridianConeRFDensity = zeros(numel(obj.enumeratedMeridianNames), numel(requestedEccentricities));
    
    % Compute variation along each of the enumerated meridians
    for meridianIndex = 1:numel(obj.enumeratedMeridianNames)
        [meridianConeRFSpacing(meridianIndex,:), meridianConeRFDensity(meridianIndex,:)] = ...
            obj.coneRFSpacingAndDensityAlongMeridian(requestedEccentricities, obj.enumeratedMeridianNames{meridianIndex}, ...
            eccUnits, densityUnits, ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio); 
    end
    
    % Do angular interpolation
    coneRFDensity2D = obj.interpolatedValuesFromMeridianValues(meridianConeRFDensity, requestedAngles);
    coneRFDensity2D = reshape(coneRFDensity2D, [eccSamplesNum eccSamplesNum]);
    midPoint = (size(coneRFDensity2D,1)-1)/2+1;
    
    
    % Transform (flip) the computed RightEye Visual Field 2D density based on requested view
    switch (theReturnedView)
        case obj.rightEyeVisualField
            % Default Watson's view, so no flipping
            xLabelString = sprintf('<-   nasal  ------- RE visual field (%s) -------  temporal ->', eccUnits);
            yLabelString = sprintf('<- inferior ------- RE visual field (%s) -------  superior ->', eccUnits);
            meridianDensities.temporal = squeeze(coneRFDensity2D(midPoint, midPoint:end));  % pos X (0 deg)
            meridianDensities.nasal    = squeeze(coneRFDensity2D(midPoint, midPoint:-1:1)); % neg X (180 deg)
            meridianDensities.superior = squeeze(coneRFDensity2D(midPoint:end, midPoint));  % pos Y (90 deg)
            meridianDensities.inferior = squeeze(coneRFDensity2D(midPoint:-1:1, midPoint)); % neg Y (270 deg)
            meridianDensities.ecc = ecc(midPoint:end);
            
        case obj.rightEyeRetina
            % Right retina view, left/right &  upside-down flip 
            coneRFDensity2D = fliplr(flipud(coneRFDensity2D));
            xLabelString = sprintf('<-   nasal  -----------  right retina (%s) -----------  temporal ->', eccUnits);
            yLabelString = sprintf('<- inferior -----------  right retina (%s) -----------  superior ->', eccUnits);
            meridianDensities.temporal = squeeze(coneRFDensity2D(midPoint, midPoint:end));  % pos X (0 deg)
            meridianDensities.nasal    = squeeze(coneRFDensity2D(midPoint, midPoint:-1:1)); % neg X (180 deg)
            meridianDensities.superior = squeeze(coneRFDensity2D(midPoint:end, midPoint));  % pos Y (90 deg)
            meridianDensities.inferior = squeeze(coneRFDensity2D(midPoint:-1:1, midPoint)); % neg Y (270 deg)
            meridianDensities.ecc = ecc(midPoint:end);
            
        case obj.leftEyeRetina
            % Left retina view, upside-down flip 
            xLabelString = sprintf('<- temporal ----------- left retina (%s) -----------  nasal ->', eccUnits);
            yLabelString = sprintf('<- inferior ----------- left retina (%s) -----------  superior ->', eccUnits);
            coneRFDensity2D = flipud(coneRFDensity2D);
            meridianDensities.temporal = squeeze(coneRFDensity2D(midPoint, midPoint:-1:1)); % neg X (180 deg)
            meridianDensities.nasal    = squeeze(coneRFDensity2D(midPoint, midPoint:end));  % pos X (0 deg)
            meridianDensities.superior = squeeze(coneRFDensity2D(midPoint:end, midPoint));  % pos Y (90 deg)
            meridianDensities.inferior = squeeze(coneRFDensity2D(midPoint:-1:1, midPoint)); % neg Y (270 deg)
            meridianDensities.ecc = ecc(midPoint:end);
    end
  
    
end
