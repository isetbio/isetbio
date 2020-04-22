function [mosaicRadiusRetinalMicrons, lambdaMicrons] = determineMosaicWidthAndLambda(mosaicWidthDegs, neuronalType)
    % Instantiate a WatsonRGCModel
    WatsonOBJ = WatsonRGCModel();
    
    % Determine smallest spacing (delta)
    switch (neuronalType)
        case 'cone'
            lambdaConesMM = WatsonOBJ.coneRFSpacingAndDensityAlongMeridian(0, 'nasal meridian', 'deg', 'mm^2', ...
                'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
            lambdaMicrons = lambdaConesMM * 1000;
        case 'mRGC'
            lambdaMidgetsMM = WatsonOBJ.mRGCRFSpacingAndDensityAlongMeridian(0, 'nasal meridian', 'deg', 'mm^2', ...
                 'adjustForISETBioConeDensity', true);
            lambdaMidgetsOnOrOffMM = sqrt(2.0)*lambdaMidgetsMM;
            lambdaMicrons = lambdaMidgetsOnOrOffMM * 1000;
        otherwise
            error('Unknown neuronalType: ''%s''.', neuronalType)
    end
    
    % Determine mosaicWidth in retinal microns
    mosaicRadiusRetinalMicrons = WatsonOBJ.rhoDegsToMMs(0.5*mosaicWidthDegs)*1000;
end

