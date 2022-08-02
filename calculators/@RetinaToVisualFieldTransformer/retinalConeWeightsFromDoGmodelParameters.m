function pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGmodelParameters(...
           retinalRFDoGparams, conesNumPooledByTheRFcenter, ...
           spatialSupportX, spatialSupportY, cm)

    % Extract fitted params
    Kc = retinalRFDoGparams.Kc;
    RcDegs = retinalRFDoGparams.RcDegs;
    RsDegs = retinalRFDoGparams.surroundToCenterRcRatio * RcDegs;
    Ks     = retinalRFDoGparams.surroundToCenterIntegratedRatio/((RsDegs/RcDegs)^2) * Kc;

    % Determine cones contributing to the RF
    coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));

    % Only consider cones whose distance from the center is < RsDegs*4
    consideredConeIndices = find(coneDistancesFromMosaicCenter < RsDegs*4);
    theConePositionsDegs = cm.coneRFpositionsDegs(consideredConeIndices,:);
    coneDistancesFromMosaicCenter = coneDistancesFromMosaicCenter(consideredConeIndices);

    [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');
    theConePositionsDegs = theConePositionsDegs(idx,:);
    consideredConeIndices = consideredConeIndices(idx);


    % Estimate RF center from the central-most conesNumPooledByTheRFcenter
    centerConeIndices = 1:conesNumPooledByTheRFcenter;
    RFcenterPos = mean(theConePositionsDegs(centerConeIndices,:),1);

    


    % Compute center cone weights
    centerConeWeights = exp(-((theConePositionsDegs(:,1)-RFcenterPos(1))/RcDegs).^2) .* exp(-((theConePositionsDegs(:,2)-RFcenterPos(2))/RcDegs).^2);
    centerConeWeights = centerConeWeights(centerConeIndices);

    % Compute surround cone weights
    minConeWeight = 0.001;
    
    surroundConeWeights = exp(-((theConePositionsDegs(:,1)-RFcenterPos(1))/RsDegs).^2) .* exp(-((theConePositionsDegs(:,2)-RFcenterPos(2))/RsDegs).^2);
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);
    
    if (~isempty(spatialSupportX))
        % Generate spatial support for computing the RFcenter/surround cone pooling maps
        [Xdegs,Ydegs] = meshgrid(spatialSupportX/60, spatialSupportY/60);
        Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

        % Compute continous center and surround weights
        centerContinuousWeights = exp(-(Rdegs/RcDegs).^2);
        surroundContinuousWeights = exp(-(Rdegs/RsDegs).^2);
        centerContinuousWeights = centerContinuousWeights(centerContinuousWeights>minConeWeight);
        surroundContinuousWeights = surroundContinuousWeights(surroundContinuousWeights>minConeWeight);

        % Compute total volume (considering cone aperture) of center & surround
        % cone weights
        retinalRFcenter2D = RetinaToVisualFieldTransformer.computeRetinalRFfromWeightsAndApertures(...
            consideredConeIndices, centerConeIndices, centerConeWeights, ...
            cm.coneApertureDiametersDegs,cm.coneRFpositionsDegs, ...
            Xdegs, Ydegs, RFcenterPos, minConeWeight);
    
        retinalRFsurround2D = RetinaToVisualFieldTransformer.computeRetinalRFfromWeightsAndApertures(...
            consideredConeIndices, surroundConeIndices, surroundConeWeights, ...
            cm.coneApertureDiametersDegs,cm.coneRFpositionsDegs, ...
            Xdegs, Ydegs, RFcenterPos, minConeWeight);

        % Match volume of center cone weights to volume of continous center weights
        correctionFactorForCenter = sum(centerContinuousWeights(:)) / sum(retinalRFcenter2D(:));
        correctionFactorForSurround = sum(surroundContinuousWeights(:)) / sum(retinalRFsurround2D(:));
    else
        correctionFactorForCenter = 1;
        correctionFactorForSurround = 1;
    end

    centerConeWeights = Kc * correctionFactorForCenter * centerConeWeights;
    surroundConeWeights = Ks * correctionFactorForSurround * surroundConeWeights;

    centerConeIndices = consideredConeIndices(centerConeIndices);
    surroundConeIndices = consideredConeIndices(surroundConeIndices);

    pooledConeIndicesAndWeightsStruct.centerConeIndices = centerConeIndices;
    pooledConeIndicesAndWeightsStruct.centerConeWeights = centerConeWeights;
    pooledConeIndicesAndWeightsStruct.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeightsStruct.surroundConeWeights = surroundConeWeights;
end