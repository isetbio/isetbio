function pooledConeIndicesAndWeights = PackerDaceyDoubleExponentModel(modelConstants, conePoolingParamsVector)

    % Retrieve params
    Kc = conePoolingParamsVector(1);
    KsToKcPeakRatio = conePoolingParamsVector(2);
    narrowToWideVolumeRatio = conePoolingParamsVector(3); % range [0.2 to 1.0]
    RwideDegs = conePoolingParamsVector(4);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = conePoolingParamsVector(5);     % range [0.01 to 1.0]

    % Compute Kwide, Knarrow, and RnarrowDegs
    [Kwide, Knarrow, RnarrowDegs] = doubleExponentSurroundRFparams(...
        Kc, RwideDegs, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);

    % This should only be executed when introducing variance to single cell surrounds, 
    % which is during the generation of the computeReadyMosaic
    % When fitting the surround, modelConstants.surroundVarianceInComputeReadyMosaic should be equal to an empty struct()
    % in which case the following block will not be executed
    if (isfield(modelConstants, 'surroundVarianceInComputeReadyMosaic')) && ...
       (~isempty(fieldnames(modelConstants.surroundVarianceInComputeReadyMosaic)))

        % Random gain with desired mean (bias) and sigma
        gaussianRandomNumber = modelConstants.surroundVarianceInComputeReadyMosaic.randomVector(1);
        gainBiasRandom  = modelConstants.surroundVarianceInComputeReadyMosaic.intSensRatioBias + ...
                           gaussianRandomNumber * modelConstants.surroundVarianceInComputeReadyMosaic.intSensRatioSigma;

        Kwide = Kwide * gainBiasRandom;
        Knarrow = Knarrow * gainBiasRandom;
    end

    % Compute center cone indices and weights
    pooledConeIndicesAndWeights.centerConeIndices = modelConstants.cachedData.centerConeIndices;
    pooledConeIndicesAndWeights.centerConeWeights = Kc * modelConstants.cachedData.centerConeWeights;

    % Compute surround cone weights from the current model params
    surroundConnectableConeWeights = ...
        Kwide * exp(-2.3*modelConstants.cachedData.distancesOfConnectableConesWithinMaxSurroundSupport/RwideDegs) + ...
        Knarrow * exp(-2.3*modelConstants.cachedData.distancesOfConnectableConesWithinMaxSurroundSupport/RnarrowDegs);
    surroundUnconnectableConeWeights = ...
        Kwide * exp(-2.3*modelConstants.cachedData.distancesOfUnconnectableConesWithinMaxSurroundSupport/RwideDegs) + ...
        Knarrow * exp(-2.3*modelConstants.cachedData.distancesOfUnconnectableConesWithinMaxSurroundSupport/RnarrowDegs);

    % Decide which surround cones to keep based on their pooling weights
    minSurroundConeWeight = mRGCMosaic.minSurroundWeightForInclusionInComputing * (Knarrow+Kwide);

    % Only keep up surround cones up to minSurroundConeWeight
    idx = find(surroundConnectableConeWeights>minSurroundConeWeight);
    pooledConeIndicesAndWeights.surroundConeIndices = ...
        modelConstants.cachedData.indicesOfConnectableConesWithinMaxSurroundSupport(idx);
    pooledConeIndicesAndWeights.surroundConeWeights = ...
        surroundConnectableConeWeights(idx);

    idx = find(surroundUnconnectableConeWeights>minSurroundConeWeight);
    pooledConeIndicesAndWeights.unconnectableSurroundConeIndices = ...
        modelConstants.cachedData.indicesOfUnconnectableConesWithinMaxSurroundSupport(idx);
    pooledConeIndicesAndWeights.unonnectableSurroundConeWeights  = ...
        surroundUnconnectableConeWeights(idx);
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
