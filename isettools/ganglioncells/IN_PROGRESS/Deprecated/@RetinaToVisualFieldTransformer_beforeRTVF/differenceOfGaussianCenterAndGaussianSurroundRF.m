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
    xDegs = modelConstants.spatialSupportDegs(:,1)-xyPos(1);
    yDegs = modelConstants.spatialSupportDegs(:,2)-xyPos(2);
    [Xdegs, Ydegs] = meshgrid(xDegs, yDegs);
    centerRF = Kc * exp(-(Xdegs/RcDegs).^2) .* exp(-(Ydegs/RcDegs).^2);
   
    % Compute the surround RF
    surroundRF = Ks * exp(-(Xdegs/RsDegs).^2) .* exp(-(Ydegs/RsDegs).^2);

    % Compute the composite RF
    theRF =  centerRF - surroundRF;
end