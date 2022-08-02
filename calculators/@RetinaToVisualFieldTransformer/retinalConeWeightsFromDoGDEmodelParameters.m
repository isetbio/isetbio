function pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGDEmodelParameters( ...
    retinalRFDoGDEparams, conesNumPooledByTheRFcenter, ...
           spatialSupportX, spatialSupportY, cm)


    % Extract fitted params
    Kc = retinalRFDoGDEparams.Kc; 
    RcDegs = retinalRFDoGDEparams.RcDegs; 
    KsToKcPeakRatio = retinalRFDoGDEparams.KsToKcPeakRatio; 
    narrowToWideVolumeRatio = retinalRFDoGDEparams.narrowToWideVolumeRatio; 
    RwideDegs = retinalRFDoGDEparams.RwideDegs; 
    RnarrowToRwideRatio = retinalRFDoGDEparams.RnarrowToRwideRatio; 

    % Compute remaining params
    RnarrowDegs = RwideDegs * RnarrowToRwideRatio;
    Kwide = (RnarrowDegs/RwideDegs) / narrowToWideVolumeRatio;
    Ks = (KsToKcPeakRatio * Kc)/(1+Kwide);

    % Determine cones contributing to the RF
    coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));

    % Only consider cones whose distance from the center is < RRwideDegs*3
    consideredConeIndices = find(coneDistancesFromMosaicCenter < RwideDegs*3);
    theConePositionsDegs = cm.coneRFpositionsDegs(consideredConeIndices,:);
    coneDistancesFromMosaicCenter = coneDistancesFromMosaicCenter(consideredConeIndices);

    [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');
    theConePositionsDegs = theConePositionsDegs(idx,:);
    consideredConeIndices = consideredConeIndices(idx);

    % Estimate RF center from the central-most conesNumPooledByTheRFcenter
    centerConeIndices = 1:conesNumPooledByTheRFcenter;
    RFcenterPos = mean(theConePositionsDegs(centerConeIndices,:),1);


    % Compute center cone weights
    centerConeWeights = exp(-((theConePositionsDegs(:,1)-RFcenterPos(1))/RcDegs).^2) .* ...
                        exp(-((theConePositionsDegs(:,2)-RFcenterPos(2))/RcDegs).^2);
    centerConeWeights = centerConeWeights(centerConeIndices);

    % Compute surround cone weights
    minConeWeight = 0.001;
    xx = abs(theConePositionsDegs(:,1) - RFcenterPos(1));
    yy = abs(theConePositionsDegs(:,2) - RFcenterPos(2));
    surroundConeWeights = exp(-2.3*xx/RnarrowDegs) .* exp(-2.3*yy/RnarrowDegs) + Kwide * ( ...
                          exp(-2.3*xx/RwideDegs) .* exp(-2.3*yy/RwideDegs));
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);
   
    if (~isempty(spatialSupportX))
        % Generate spatial support for computing the RFcenter/surround cone pooling maps
        [Xdegs,Ydegs] = meshgrid(spatialSupportX/60, spatialSupportY/60);
        Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

        % Compute continous center and surround weights
        centerContinuousWeights = exp(-(Rdegs/RcDegs).^2);
        surroundContinuousWeights = exp(-2.3*Rdegs/RnarrowDegs) + Kwide * exp(-2.3*Rdegs/RwideDegs);
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