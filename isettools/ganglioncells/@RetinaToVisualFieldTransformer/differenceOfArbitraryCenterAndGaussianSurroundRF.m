function theRF = differenceOfArbitraryCenterAndGaussianSurroundRF(...
       RFcenterConeMap, RFcenterCharacteristicRadiusDegs, paramsVector, spatialSupportDegs)

    % Extract params for the fit
    Kc = paramsVector(1);
    surroundToCenterCharacteristicRadiiRatio = paramsVector(2);
    surroundToCenterIntegratedSensitivitiesRatio = paramsVector(3);
    RsDegs = surroundToCenterCharacteristicRadiiRatio * RFcenterCharacteristicRadiusDegs;
    Ks = Kc * surroundToCenterIntegratedSensitivitiesRatio/(surroundToCenterCharacteristicRadiiRatio^2);

    % Spatial support mesh
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the center RF
    centerRF = Kc * RFcenterConeMap;

    % Compute the surround RF
    surroundRF = Ks * exp(-Rdegs2/RsDegs^2);

    % Correct so as to obtain the desired  surroundToCenterIntegratedSensitivitiesRatio
    surroundBoost = surroundToCenterIntegratedSensitivitiesRatio/(sum(surroundRF(:))/sum(centerRF(:)));
    surroundRF = surroundRF * surroundBoost;

    % Compute the composite RF
    theRF = centerRF - surroundRF;

    %actualSCintegratedSensitivityRatio = sum(surroundRF(:))/sum(centerRF(:));
end