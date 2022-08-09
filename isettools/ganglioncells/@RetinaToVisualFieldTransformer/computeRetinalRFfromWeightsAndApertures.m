function retinalRF = computeRetinalRFfromWeightsAndApertures(...
        consideredConeIndices, coneIndices, coneWeights, ...
        coneApertureDiametersDegs, coneRFpositionsDegs, ...
        Xdegs, Ydegs, RFcenter, minConeWeight)

    
    Xdegs = Xdegs + RFcenter(1);
    Ydegs = Ydegs + RFcenter(2);

    for iCone = 1:numel(coneIndices)
        theConeIndex = consideredConeIndices(coneIndices(iCone));
        coneRcDegs = 0.204 * sqrt(2.0) * coneApertureDiametersDegs(theConeIndex);

        % Gaussian-shaped aperture
        xx = (Xdegs-coneRFpositionsDegs(theConeIndex,1))/coneRcDegs;
        yy = (Ydegs-coneRFpositionsDegs(theConeIndex,2))/coneRcDegs;
        theConeApertureRF = exp(-(xx.^2 + yy.^2));
        theConeApertureRF(theConeApertureRF<minConeWeight) = 0;

        if (iCone == 1)
            retinalRF = coneWeights(iCone) * theConeApertureRF;
        else
            retinalRF = retinalRF + coneWeights(iCone) * theConeApertureRF;
        end
    end
end