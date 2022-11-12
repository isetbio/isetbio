function [theRF, centerRF, surroundRF] = differenceOfEllipsoidalGaussianCenterAndGaussianSurroundRF(...
    modelConstants, paramsVector)

    % Extract params for the fit
    Kc = paramsVector(1);
    RcDegsX = paramsVector(2);
    RcDegsY = paramsVector(3);
    xyPos(1) = real(paramsVector(4));
    xyPos(2) = imag(paramsVector(4));
    exponentX = paramsVector(5);
    exponentY = paramsVector(6);
    rotationDegs = paramsVector(7);
    surroundToCenterCharacteristicRadiiRatio = paramsVector(8);
    surroundToCenterIntegratedSensitivitiesRatio = paramsVector(9);
    RsDegs = surroundToCenterCharacteristicRadiiRatio * sqrt(RcDegsX*RcDegsY);
    Ks = Kc * surroundToCenterIntegratedSensitivitiesRatio/(surroundToCenterCharacteristicRadiiRatio^2);

    % Compute the center RF
    xDegs = modelConstants.spatialSupportDegs(:,1)-xyPos(1);
    yDegs = modelConstants.spatialSupportDegs(:,2)-xyPos(2);
    [Xdegs, Ydegs] = meshgrid(xDegs, yDegs);
    centerRF = Kc * exp(-(abs(Xdegs)/RcDegsX).^(2*exponentX)) .* exp(-(abs(Ydegs)/RcDegsY).^(2*exponentY));
    centerRF = imrotate(centerRF, rotationDegs, "bilinear", "crop");

    % Compute the surround RF
    surroundRF = Ks * exp(-(Xdegs/RsDegs).^2) .* exp(-(Ydegs/RsDegs).^2);

    % Compute the composite RF
    theRF =  centerRF - surroundRF;
end