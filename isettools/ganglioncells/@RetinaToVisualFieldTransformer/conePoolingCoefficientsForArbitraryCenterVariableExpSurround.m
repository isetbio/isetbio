function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterVariableExpSurround(...
            modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    Ks = conePoolingParamsVector(2);
    % Params for the variable exponentials
    RsDegs = conePoolingParamsVector(3);               
    variableExponent = conePoolingParamsVector(4); 
    
    % compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2));

    % Compute surround cone weights
    minConeWeight = 0.001;

    surroundConeWeights = exp(-(abs(coneDistancesFromRFCenter/RsDegs)).^variableExponent);
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = Ks * surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
end
