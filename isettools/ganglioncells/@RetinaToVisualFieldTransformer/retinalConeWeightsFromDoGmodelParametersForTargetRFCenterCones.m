function pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGmodelParametersForTargetRFCenterCones(...
           retinalRFDoGparams, cm, targetRFCenterConesIndices)


    % Extract fitted params
    Kc = retinalRFDoGparams.Kc;
    RcDegs = retinalRFDoGparams.RcDegs;
    RsDegs = retinalRFDoGparams.surroundToCenterRcRatio * RcDegs;
    Ks     = retinalRFDoGparams.surroundToCenterIntegratedRatio/((RsDegs/RcDegs)^2) * Kc;

    % Estimate RF center from the targetRFCenterConesIndices
    centerConePositions = cm.coneRFpositionsDegs(targetRFCenterConesIndices,:);
    RFcenterPos = mean(centerConePositions,1);

    % Compute center cone weights
    centerConeWeights = exp(-((centerConePositions(:,1)-RFcenterPos(1))/RcDegs).^2) .* ...
                        exp(-((centerConePositions(:,2)-RFcenterPos(2))/RcDegs).^2);

    % Compute surround cone weights
    surroundConeWeights = exp(-((cm.coneRFpositionsDegs(:,1)-RFcenterPos(1))/RsDegs).^2) .* ...
                          exp(-((cm.coneRFpositionsDegs(:,2)-RFcenterPos(2))/RsDegs).^2);
    minConeWeight = 0.001;
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeightsStruct.centerConeIndices = targetRFCenterConesIndices;
    pooledConeIndicesAndWeightsStruct.surroundConeIndices = surroundConeIndices;

    pooledConeIndicesAndWeightsStruct.centerConeWeights = Kc*centerConeWeights;
    pooledConeIndicesAndWeightsStruct.surroundConeWeights = Ks*surroundConeWeights;
end