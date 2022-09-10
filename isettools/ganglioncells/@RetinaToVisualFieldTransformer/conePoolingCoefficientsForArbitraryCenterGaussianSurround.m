function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianSurround(...
    modelConstants, conePoolingParamsVector)

    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = modelConstants.weightsOfCenterCones;

    % compute surround cone indices and weights
    Ks = conePoolingParamsVector(1);
    RsDegs = conePoolingParamsVector(2);

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);

    % Compute surround cone weights
    minConeWeight = 0.005;
    
    surroundConeWeights = exp(-coneDistancesFromRFCenterSquared/RsDegs.^2);
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = Ks * surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
end