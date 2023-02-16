function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround( ...
    modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]

    % Compute Kwide, Knarrow, and RnarrowDegs
    [Kwide, Knarrow, RnarrowDegs] = RTVF.H1doubleExponentRFparams(...
        Kc, RwideDegs, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);


    % Compute center cone indices and weights
    centerConeIndices = modelConstants.indicesOfCenterCones;
    centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center position
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(centerConeIndices,:),1);

    % Compute the distances of ALL cones in the input cone mosaic from the RF center - THIS SHOULD BE CACHED
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);
    coneDistancesFromRFCenter = sqrt(coneDistancesFromRFCenterSquared);

    % Note: some of these surround cones will be outside of the RF spatial support, so
    % eliminated them here
    surroundConeIndices = find(coneDistancesFromRFCenterSquared <= max(modelConstants.Rdegs2(:)));
    coneDistancesFromRFCenter = coneDistancesFromRFCenter(surroundConeIndices);

    % Compute surround cone weights from the current model params
    surroundConeWeights = Kwide * exp(-2.3*coneDistancesFromRFCenter/RwideDegs) + ...
                          Knarrow * exp(-2.3*coneDistancesFromRFCenter/RnarrowDegs);

    % Decide which surround cones to keep based on their pooling weights
    minSurroundConeWeight = 1e-3 * (Knarrow+Kwide);
    idx = find(surroundConeWeights>minSurroundConeWeight);
    surroundConeIndices = surroundConeIndices(idx);
    surroundConeWeights = surroundConeWeights(idx);

    % Decide which surround cones to keep based on their type
    [surroundConeIndices, surroundConeWeights, ...
     nonConnectableSurroundConeIndices, ...
     nonConnectableSurroundConeWeights] = RTVF.connectableSurroundConeIndicesAndWeights(...
         surroundConeIndices, surroundConeWeights, modelConstants);

    % The indices of center and surround cones
    pooledConeIndicesAndWeights.centerConeIndices = centerConeIndices;
    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.nonConnectableSurroundConeIndices = nonConnectableSurroundConeIndices;

    % Unlike the old version, we do not need compensate for variations in relative cone efficiency
    % resulting from c variations in cone IS diameter and cone OS length with eccentricity 
    % because we now take this into account in the cone weights to compensate for variations in
    % the RTV.retinalSubregionConeMapFromPooledConeInputs() method
    pooledConeIndicesAndWeights.centerConeWeights = centerConeWeights;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
    pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights = nonConnectableSurroundConeWeights;
end    