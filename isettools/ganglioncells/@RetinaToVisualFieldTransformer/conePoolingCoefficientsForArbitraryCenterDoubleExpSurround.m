function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround(...
            modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    % Params for the sum of 2 exponentials (wide field + narrow field)
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]
    
    RnarrowDegs = RwideDegs * RnarrowToRwideRatio;
    Kwide = (RnarrowDegs/RwideDegs) / narrowToWideVolumeRatio;
    Ks = (KsToKcPeakRatio * Kc)/(1+Kwide);

    % compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2));

    % Compute surround cone weights
    minConeWeight = 0.001*(1+Kwide);

    surroundConeWeights = exp(-2.3*coneDistancesFromRFCenter/RnarrowDegs) + Kwide * exp(-2.3*coneDistancesFromRFCenter/RwideDegs);
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = Ks * surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
end
