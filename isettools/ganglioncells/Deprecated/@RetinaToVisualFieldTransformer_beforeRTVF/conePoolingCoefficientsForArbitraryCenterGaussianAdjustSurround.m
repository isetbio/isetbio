function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianAdjustSurround(...
    modelConstants, conePoolingParamsVector)
    
    % Retrieve params
    Kc = conePoolingParamsVector(1);
    Ks = conePoolingParamsVector(2);
    RsDegs = conePoolingParamsVector(3);
    surroundRadialAdjustmentsVector = conePoolingParamsVector(4:end);

    % compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);

    % Compute surround cone weights
    minConeWeight = 0.001;
    
    % adjustment for some surround cone weights
    surroundConeIndicesToReceiveWeightsAdjustments = ...
        find(coneDistancesFromRFCenterSquared < (max(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs))^2);

    surroundConeRadiiToReceiveWeightsAdjustments = sqrt(coneDistancesFromRFCenterSquared(surroundConeIndicesToReceiveWeightsAdjustments));
    surroundAdjustmentFactors = interp1(...
        modelConstants.arbitrarySurroundCorrectionRadialSupportDegs, surroundRadialAdjustmentsVector, surroundConeRadiiToReceiveWeightsAdjustments, ...
        'linear', 'extrap');

    % Gaussian surround weights
    surroundConeWeights = exp(-coneDistancesFromRFCenterSquared/(RsDegs.^2));

    % Adjust weights in certain surround cones
    surroundConeWeights(surroundConeIndicesToReceiveWeightsAdjustments) = ...
        surroundConeWeights(surroundConeIndicesToReceiveWeightsAdjustments) + surroundAdjustmentFactors/Ks;

    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = Ks * surroundConeWeights(surroundConeIndices);
   
    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;

    error('This function needs to be updated by calling RetinaToVisualFieldTransformer.connectableSurroundConeIndicesAndWeights and RetinaToVisualFieldTransformer.coneEfficacyAdjustedGains')
end