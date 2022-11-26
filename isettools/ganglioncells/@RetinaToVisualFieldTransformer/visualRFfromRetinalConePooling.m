function [theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
    visualRFfromRetinalConePooling(modelConstants, retinalPoolingParams)

    % Compute the indices and retinal weights of cones pooled by the center and the
    % surround mechanism
    pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(modelConstants, retinalPoolingParams);

    % Compute theRetinalRFcenterConeMap
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.centerConeIndices,:);
    rfCenterPositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        modelConstants.theConeMosaic,...
        theConePositionsDegs, ...
        pooledConeIndicesAndWeights.centerConeWeights, ...
        modelConstants.spatialSupportDegs);

    % Compute theRetinalRFsurroundConeMap
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.surroundConeIndices,:);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        modelConstants.theConeMosaic,...
        theConePositionsDegs, ...
        pooledConeIndicesAndWeights.surroundConeWeights, ...
        modelConstants.spatialSupportDegs);


    % The composite retinal RF cone map
    theRetinalRF = theRetinalRFcenterConeMap - theRetinalRFsurroundConeMap;

    % Convolve the composite cone-pooling based retinal RF to get the corresponding visual RF
    theVisualRF = conv2(theRetinalRF, modelConstants.thePSF, 'same');
end