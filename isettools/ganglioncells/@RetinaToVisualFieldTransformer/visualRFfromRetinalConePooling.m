function [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = visualRFfromRetinalConePooling(modelConstants, retinalPoolingParams)

    % Compute the indices and retinal weights of cones pooled by the center and the
    % surround mechanism
    pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(modelConstants, retinalPoolingParams);

    % Compute theRetinalRFcenterConeMap
    theConeApertureDiameterDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(pooledConeIndicesAndWeights.centerConeIndices);
    theConeCharacteristicRadiiDegs = ...
        modelConstants.coneCharacteristicRadiusConversionFactor * theConeApertureDiameterDegs;
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.centerConeIndices,:);
    rfCenterPositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConeApertureDiameterDegs, theConePositionsDegs, pooledConeIndicesAndWeights.centerConeWeights, ...
        modelConstants.spatialSupportDegs);


    % Compute theRetinalRFsurroundConeMap
    theConeApertureDiameterDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(pooledConeIndicesAndWeights.surroundConeIndices);
    theConeCharacteristicRadiiDegs = ...
        modelConstants.coneCharacteristicRadiusConversionFactor * theConeApertureDiameterDegs;
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.surroundConeIndices,:);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConeApertureDiameterDegs, theConePositionsDegs, pooledConeIndicesAndWeights.surroundConeWeights, ...
        modelConstants.spatialSupportDegs);


    % The composite retinal RF cone map
    theRetinalRF = theRetinalRFcenterConeMap - theRetinalRFsurroundConeMap;

    % Convolve the composite cone-pooling based retinal RF to get the corresponding visual RF
    theFittedVisualRF = conv2(theRetinalRF, modelConstants.thePSF, 'same');

end