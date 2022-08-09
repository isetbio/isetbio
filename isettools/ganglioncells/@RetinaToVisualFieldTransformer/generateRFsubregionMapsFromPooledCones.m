function [retinalRFcenter2D, retinalRFsurround2D] = generateRFsubregionMapsFromPooledCones(...
        spatialSupportX, spatialSupportY, cm, pooledConeIndicesAndWeightsStruct)

    centerConeIndices = pooledConeIndicesAndWeightsStruct.centerConeIndices;
    centerConeWeights = pooledConeIndicesAndWeightsStruct.centerConeWeights;
    surroundConeIndices = pooledConeIndicesAndWeightsStruct.surroundConeIndices;
    surroundConeWeights = pooledConeIndicesAndWeightsStruct.surroundConeWeights;

    RFcenter = mean(cm.coneRFpositionsDegs(centerConeIndices,:),1);
    minConeWeight = 0.001;

    % Generate spatial support for computing the RFcenter/surround cone pooling maps
    [Xdegs,Ydegs] = meshgrid(spatialSupportX+RFcenter(1), spatialSupportY+RFcenter(2));
   
    for iCone = 1:numel(centerConeIndices)
        theConeIndex = centerConeIndices(iCone);
        coneRcDegs = 0.204 * sqrt(2.0) * cm.coneApertureDiametersDegs(theConeIndex);

        % Gaussian-shaped aperture
        xx = (Xdegs-cm.coneRFpositionsDegs(theConeIndex,1))/coneRcDegs;
        yy = (Ydegs-cm.coneRFpositionsDegs(theConeIndex,2))/coneRcDegs;
        theConeApertureRF = exp(-(xx.^2 + yy.^2));
        theConeApertureRF(theConeApertureRF<minConeWeight) = 0;

        if (iCone == 1)
            retinalRFcenter2D = theConeApertureRF * centerConeWeights(iCone);
        else
            retinalRFcenter2D = retinalRFcenter2D + theConeApertureRF* centerConeWeights(iCone);
        end
    end


    for iCone = 1:numel(surroundConeIndices)
        theConeIndex = surroundConeIndices(iCone);
        coneRcDegs = 0.204 * sqrt(2.0) * cm.coneApertureDiametersDegs(theConeIndex);

        % Gaussian-shaped aperture
        xx = (Xdegs-cm.coneRFpositionsDegs(theConeIndex,1))/coneRcDegs;
        yy = (Ydegs-cm.coneRFpositionsDegs(theConeIndex,2))/coneRcDegs;
        theConeApertureRF = exp(-(xx.^2 + yy.^2));
        theConeApertureRF(theConeApertureRF<minConeWeight) = 0;

        if (iCone == 1)
            retinalRFsurround2D = theConeApertureRF * surroundConeWeights(iCone);
        else
            retinalRFsurround2D = retinalRFsurround2D + theConeApertureRF * surroundConeWeights(iCone); 
        end
    end 
end