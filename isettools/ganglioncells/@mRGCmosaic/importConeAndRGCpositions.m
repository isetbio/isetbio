function [coneRFpositionsMicrons, coneRFpositionsDegs, ...
          rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = importConeAndRGCpositions(eccentricityDegs, sizeDegs, whichEye)

    % Import mRGC RF positions for passed eccentricity, size, and eye
    [rgcRFpositionsMicrons,  rgcRFpositionsDegs] = retinalattice.import.mRGCPositions(eccentricityDegs, sizeDegs, whichEye);
    
    % Determine mean surround subregion size (2 sigma) for this eccentricity
    maxEccentricityDegs = mean(sqrt(sum(rgcRFpositionsDegs.^2,2)));
    extraDegsForRGCSurround = 2.0 * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);
    
    % Import cone RF positions, for passed eccentricity, size, and eye
    [coneRFpositionsMicrons, coneRFpositionsDegs] = retinalattice.import.conePositions(eccentricityDegs, ...
        sizeDegs, whichEye);
                           
end


