function val = totalRGCRFDensityAlongMeridian(obj, eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits)

    % Validate eccUnits
    obj.validateEccUnits(eccUnits);
    
    % Validate densityUnits
    obj.validateDensityUnits(densityUnits);
    
    % Make sure eccentricities is a 1xN vector
    if (size(eccentricities,1)>1)
        eccentricities = eccentricities';
    end
    assert(size(eccentricities,1) == 1, 'Eccentricities must be a 1xN vector');
    
    % Convert passed eccentricities as requested
    switch (eccUnits)
        case obj.visualDegsEccUnits
            eccDegs = eccentricities;
        case obj.retinalMMEccUnits
            eccDegs = obj.rhoMMsToDegs(eccentricities);
    end
    
    % Step 1. compute peak total RGC RF density 
    % Retrieve peak cone density
    [~,peakConeDensity] = obj.coneRFSpacingAndDensityAlongMeridian(0,'temporal meridian', 'deg', densityUnits);


    % The foveal density of midget RGC RFs is twice the cone density
    % because in the fovea each cone connects to exactly 2 midget RGCs (one
    % ON and one OFF). This is Equation (1) in the Watson (2014) paper.
    peakMidgetRGCRFDensity = 2 * peakConeDensity;
    
    % Retrieve percentage of total ganglion cells that are midget at 0 eccentricity
    percentageOfRGCsThatAreMidgetsAtZeroEccentricity = obj.f0;
    % This is Equation (2) in the Watson (2014) paper.
    peakTotalRGCRFDensity = 1/percentageOfRGCsThatAreMidgetsAtZeroEccentricity * peakMidgetRGCRFDensity;
    
    % Step 2. Get the meridian params for the requested meridian
    meridianParams = obj.meridianParams(rightEyeVisualFieldMeridianName);
    
    % This is equation (4) in the Watson (2014) paper.
    eccVariation = ...
        meridianParams.a_k * (1+eccDegs/meridianParams.r_2k).^(-2) + ...
        (1-meridianParams.a_k) * exp(-eccDegs/meridianParams.r_ek);

    % Return the total RGC RF density along the requested meridian
    val = peakTotalRGCRFDensity * eccVariation;
end

