function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround(...
            modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]


    % Compute Kwide, Knarrow, and RnarrowDegs
    [Kwide, Knarrow, RnarrowDegs] = RetinaToVisualFieldTransformer.H1doubleExponentRFparams(...
        Kc, RwideDegs, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);


    % Compute center cone indices and weights
    centerConeIndices = modelConstants.indicesOfCenterCones;
    centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    
    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2));

    % Compute surround cone weights from the current model params
    surroundConeWeights = Kwide * exp(-2.3*coneDistancesFromRFCenter/RwideDegs) + Knarrow * exp(-2.3*coneDistancesFromRFCenter/RnarrowDegs);


    minSurroundConeWeight = 1e-3 * (Knarrow+Kwide);
    surroundConeIndices = find(surroundConeWeights>minSurroundConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);

    % Keep only the connectable surround cones
    [surroundConeIndices, surroundConeWeights, ...
     nonConnectableSurroundConeIndices, ...
     nonConnectableSurroundConeWeights] = RetinaToVisualFieldTransformer.connectableSurroundConeIndicesAndWeights(...
         surroundConeIndices, surroundConeWeights, modelConstants);

    % The indices of center and surround cones
    pooledConeIndicesAndWeights.centerConeIndices = centerConeIndices;
    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.nonConnectableSurroundConeIndices = nonConnectableSurroundConeIndices;
    
    if (modelConstants.coneWeightsCompensateForVariationsInConeEfficiency)
        % Adjust cone weights to compensate for variations in
        % relative efficiency of the input cones 
        % (which result due to variations in cone IS diameter and cone OS length with eccentricity)
        maxEfficiency = [];
        [pooledConeIndicesAndWeights.centerConeWeights, maxEfficiency] = RetinaToVisualFieldTransformer.coneEfficacyAdjustedGains(...
            modelConstants.theConeMosaic, ...
            modelConstants.theConeMosaic.coneApertureDiametersDegs(centerConeIndices), ...
            modelConstants.theConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(centerConeIndices), ...
            centerConeWeights, ...
            maxEfficiency);
    
        % Adjust the surround cone weights to counteract the cone efficacy gain
        pooledConeIndicesAndWeights.surroundConeWeights = RetinaToVisualFieldTransformer.coneEfficacyAdjustedGains(...
            modelConstants.theConeMosaic, ...
            modelConstants.theConeMosaic.coneApertureDiametersDegs(surroundConeIndices), ...
            modelConstants.theConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(surroundConeIndices), ...
            surroundConeWeights, ...
            maxEfficiency);

        pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights = RetinaToVisualFieldTransformer.coneEfficacyAdjustedGains(...
            modelConstants.theConeMosaic, ...
            modelConstants.theConeMosaic.coneApertureDiametersDegs(nonConnectableSurroundConeIndices), ...
            modelConstants.theConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(nonConnectableSurroundConeIndices), ...
            nonConnectableSurroundConeWeights, ...
            maxEfficiency);
    
    else
        pooledConeIndicesAndWeights.centerConeWeights = centerConeWeights;
        pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
        pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights = nonConnectableSurroundConeWeights;
    end
end



