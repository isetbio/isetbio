function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleGaussianSurround(...
    modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    Ks = conePoolingParamsVector(2);
    RsDegs = conePoolingParamsVector(3);
    Kwide = conePoolingParamsVector(4);
    RsWideToRsRatio = conePoolingParamsVector(5);
    RsWideDegs = RsWideToRsRatio * RsDegs;

    % compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);

    % Compute surround cone weights
    minConeWeight = 0.001*(min([1 Kwide]));
    
    surroundConeWeights = exp(-coneDistancesFromRFCenterSquared/(RsDegs.^2)) + Kwide * exp(-coneDistancesFromRFCenterSquared/(RsWideDegs.^2));
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = Ks * surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
end