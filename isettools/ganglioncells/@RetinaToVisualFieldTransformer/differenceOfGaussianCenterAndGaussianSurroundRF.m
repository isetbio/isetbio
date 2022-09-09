function theRF = differenceOfGaussianCenterAndGaussianSurroundRF(paramsVector, spatialSupportDegs)

    % Extract params for the fit
    Kc = paramsVector(1);
    RcDegs = paramsVector(2);
    surroundToCenterCharacteristicRadiiRatio = paramsVector(3);
    surroundToCenterIntegratedSensitivitiesRatio = paramsVector(4);
    RsDegs = surroundToCenterCharacteristicRadiiRatio * RcDegs;
    Ks = Kc * surroundToCenterIntegratedSensitivitiesRatio/(surroundToCenterCharacteristicRadiiRatio^2);

    % Spatial support mesh
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the center RF
    centerRF = Kc * exp(-Rdegs2/RcDegs^2);

    % Compute the surround RF
    surroundRF = Ks * exp(-Rdegs2/RsDegs^2);

    % Compute the composite RF
    theRF =  centerRF - surroundRF;

    %actualSCintegratedSensitivityRatio = sum(surroundRF(:))/sum(centerRF(:))
end