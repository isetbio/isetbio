function angularEccDegs = linearEccMMsToAngularEccDegs(obj,linearEccMM)

    angularEccDegs = obj.inputConeMosaic.distanceMicronsToDistanceDegreesForCmosaic(linearEccMM * 1e3);

end