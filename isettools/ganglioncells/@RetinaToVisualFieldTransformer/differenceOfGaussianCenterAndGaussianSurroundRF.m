function [theRF, centerRF, surroundRF] = differenceOfGaussianCenterAndGaussianSurroundRF(...
    modelConstants, paramsVector)

    % Extract params for the fit
    Kc = paramsVector(1);
    RcDegs = paramsVector(2);

    % offset not taken into account yet
    xyPos(1) = real(paramsVector(3));
    xyPos(2) = imag(paramsVector(3));

    surroundToCenterCharacteristicRadiiRatio = paramsVector(4);
    surroundToCenterIntegratedSensitivitiesRatio = paramsVector(5);
    RsDegs = surroundToCenterCharacteristicRadiiRatio * RcDegs;
    Ks = Kc * surroundToCenterIntegratedSensitivitiesRatio/(surroundToCenterCharacteristicRadiiRatio^2);

    % Compute the center RF
    centerRF = Kc * exp(-modelConstants.Rdegs2/RcDegs^2);

    % Compute the surround RF
    surroundRF = Ks * exp(-modelConstants.Rdegs2/RsDegs^2);

    % Compute the composite RF
    theRF =  centerRF - surroundRF;
end