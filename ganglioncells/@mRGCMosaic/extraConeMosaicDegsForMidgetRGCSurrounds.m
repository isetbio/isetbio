function extraDegs = extraConeMosaicDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs)
    
    maxEccentricityDegs = sqrt(sum((eccentricityDegs + 0.5*sizeDegs).^2,2));
    maxSurroundSigmaDegs = 1/sqrt(2.0) * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);
    % At 2.145 x sigma, a 2D Gaussian is at 1% of the peak
    % At 2.628 x sigma, a 2D Gaussian is at 0.1% of the peak
    extraDegs = 2.146 * maxSurroundSigmaDegs;

end