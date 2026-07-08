function linearEccMMs = angularEccDegsToLinearEccMMs(obj, angularEccDegs)

    linearEccMMs = 1e-3 * obj.inputConeMosaic.distanceDegreesToDistanceMicronsForCmosaic(angularEccDegs);

end