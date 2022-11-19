function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround(...
            modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]

    % Params for the sum of 2 exponentials (wide field + narrow field, Packer& Dacey 2002)
    % RF = Kwide * exp(-2.3*R/Rwide) + Knarrow * exp(-2.3*R/Rnarrow) 
    % Rwide: radius at which sensitivity drops to 10%, which is defined as half the RFdiameter
    % Volume: R^2 * K
    % N/W volume ratio (range: 0.2 to 1.0) = Knarrow/Kwide * (Rnarrow/Rwide)^2;

    % Knarrow = volumeRatio * Kwide / (RnarrowToRwideRatio^2)
    % Rnarrow = Rwide * RnarrowToRwideRatio;

    % Kwide + Knarrow = Ks = KsToKcPeakRatio * Kc;
    % Kwide = KsToKcPeakRatio * Kc - Knarrow;
    % Kwide = KsToKcPeakRatio * Kc - volumeRatio * Kwide / (RnarrowToRwideRatio^2);
    % Kwide*(1+volumeRatio/(RnarrowToRwideRatio^2)) = KsToKcPeakRatio * Kc;
    % Kwide = KsToKcPeakRatio * Kc / (1+volumeRatio/(RnarrowToRwideRatio^2));

    % Compute parameters
    Kwide = KsToKcPeakRatio * Kc / (1+narrowToWideVolumeRatio/(RnarrowToRwideRatio^2));
    Knarrow = Kwide * narrowToWideVolumeRatio / (RnarrowToRwideRatio^2);
    RnarrowDegs = RwideDegs * RnarrowToRwideRatio;

    % Compute center cone indices and weights
    centerConeIndices = modelConstants.indicesOfCenterCones;
    centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    
    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2));

    % Compute surround cone weights
    surroundConeWeights = Kwide * exp(-2.3*coneDistancesFromRFCenter/RwideDegs) + Knarrow * exp(-2.3*coneDistancesFromRFCenter/RnarrowDegs);
    minSurroundConeWeight = 1e-3 * (Knarrow+Kwide);
    surroundConeIndices = find(surroundConeWeights>minSurroundConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);

    % Keep only the connectable surround cones
    [surroundConeIndices, surroundConeWeights] = RetinaToVisualFieldTransformer.connectableSurroundConeIndicesAndWeights(...
         surroundConeIndices, surroundConeWeights, modelConstants);

    % The indices of center and surround cones
    pooledConeIndicesAndWeights.centerConeIndices = centerConeIndices;
    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;

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
    else
        pooledConeIndicesAndWeights.centerConeWeights = centerConeWeights;
        pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
    end
end



