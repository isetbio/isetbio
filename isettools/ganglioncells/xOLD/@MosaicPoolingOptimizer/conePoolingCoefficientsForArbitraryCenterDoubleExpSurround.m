function pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround( ...
    modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]

    % Compute Kwide, Knarrow, and RnarrowDegs
    [Kwide, Knarrow, RnarrowDegs] = doubleExponentSurroundRFparams(...
        Kc, RwideDegs, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);

    % Compute center cone indices and weights
    centerConeIndices = modelConstants.indicesOfCenterCones;
    centerConeWeights = Kc * modelConstants.weightsOfCenterCones;

    surroundConeIndices = modelConstants.cachedData.surroundConeIndices;
    coneDistancesFromRFCenter = modelConstants.cachedData.coneDistancesFromRFCenter;


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
     nonConnectableSurroundConeWeights] = MosaicPoolingOptimizer.connectableSurroundConeIndicesAndWeights(...
         surroundConeIndices, surroundConeWeights, modelConstants);

    % The indices of center and surround cones
    pooledConeIndicesAndWeights.centerConeIndices = centerConeIndices;
    pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
    pooledConeIndicesAndWeights.nonConnectableSurroundConeIndices = nonConnectableSurroundConeIndices;

    pooledConeIndicesAndWeights.centerConeWeights = centerConeWeights;
    pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;
    pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights = nonConnectableSurroundConeWeights;
end   

function [Kwide, Knarrow, Rnarrow] = doubleExponentSurroundRFparams(...
    Kc, Rwide, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio)
% Given Kc (center gain), KsToKcRatio (surround/center peak sensitivity
% ratio) and Rwide (wide RF surround radius), narrowToWideVolumeRatio and RnarrowToRwideRatio,
% compute the peak sensitivity gains for the wide and the narrow surround components

    % Params for the sum of 2 exponentials (wide field + narrow field, Packer& Dacey 2002)
    % RF = Kwide * exp(-2.3*R/Rwide) + Knarrow * exp(-2.3*R/Rnarrow) 
    % Rwide: radius at which sensitivity drops to 10%, which is defined as half the RFdiameter
    % Volume: R^2 * K

    Kwide = KsToKcPeakRatio*Kc / (1 + narrowToWideVolumeRatio / (RnarrowToRwideRatio^2));
    Knarrow  = Kwide * narrowToWideVolumeRatio / (RnarrowToRwideRatio^2);
    Rnarrow  = RnarrowToRwideRatio * Rwide;
end
