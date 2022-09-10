function theRF = differenceOfGaussianCenterAndGaussianSurroundRF(...
    modelConstants, paramsVector)

    % Extract params for the fit
    Kc = paramsVector(1);
    RcDegs = paramsVector(2);
    surroundToCenterCharacteristicRadiiRatio = paramsVector(3);
    surroundToCenterIntegratedSensitivitiesRatio = paramsVector(4);
    RsDegs = surroundToCenterCharacteristicRadiiRatio * RcDegs;
    Ks = Kc * surroundToCenterIntegratedSensitivitiesRatio/(surroundToCenterCharacteristicRadiiRatio^2);

    % Compute the center RF
    centerRF = Kc * exp(-modelConstants.Rdegs2/RcDegs^2);

    % Compute the surround RF
    surroundRF = Ks * exp(-modelConstants.Rdegs2/RsDegs^2);

    % Compute the composite RF
    theRF =  centerRF - surroundRF;

    %actualSCintegratedSensitivityRatio = sum(surroundRF(:))/sum(centerRF(:))
end