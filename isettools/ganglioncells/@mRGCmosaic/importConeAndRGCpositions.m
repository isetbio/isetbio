function [coneRFpositionsMicrons, coneRFpositionsDegs, ...
          rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = importConeAndRGCpositions(sourceLatticeSizeDegs, eccentricityDegs, sizeDegs, whichEye)

    % Import mRGC RF positions for passed eccentricity, size, and eye
    [rgcRFpositionsMicrons,  rgcRFpositionsDegs] = retinalattice.import.finalMRGCPositions(sourceLatticeSizeDegs, eccentricityDegs, sizeDegs, whichEye);

    % Determine mean surround subregion size (2 sigma) for this eccentricity
    meanEccentricityDegs = mean(sqrt(sum(rgcRFpositionsDegs.^2,2)));
    extraDegsForRGCSurround = 2.0 * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(meanEccentricityDegs);
    
    % Import cone RF positions, for passed eccentricity, size, and eye
    [coneRFpositionsMicrons, coneRFpositionsDegs] = retinalattice.import.finalConePositions(sourceLatticeSizeDegs, eccentricityDegs, ...
        sizeDegs, whichEye);
                           
end


