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
    Knarrow = narrowToWideVolumeRatio * Kwide / (RnarrowToRwideRatio^2);
    RnarrowDegs = RwideDegs * RnarrowToRwideRatio;

    % compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.indicesOfCenterCones;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    % Compute the RF center
    RFcenterPos = mean(modelConstants.theConeMosaic.coneRFpositionsDegs(modelConstants.indicesOfCenterCones,:),1);

    % Compute all cone distances from the RF center (squared)
    coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, modelConstants.theConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2));

    % Compute surround cone weights
    minConeWeight = 0.01*(min([Knarrow Kwide]));

    surroundConeWeights = Kwide * exp(-2.3*coneDistancesFromRFCenter/RwideDegs) + Knarrow * exp(-2.3*coneDistancesFromRFCenter/RnarrowDegs);
    surroundConeIndices = find(surroundConeWeights>minConeWeight);
    surroundConeWeights = surroundConeWeights(surroundConeIndices);

    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
end
