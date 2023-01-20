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

    % Convolve the retinal RF center with the center cone PSF to get the corresponding visual RF center
    theVisualRFcenter = conv2(theRetinalRFcenterConeMap, modelConstants.theRFCenterConeMajorityPSF, 'same');

    % Convolve the retinal RF surround with the surround L+M PSF to get the corresponding visual RF surround
    theVisualRFsurround = conv2(theRetinalRFsurroundConeMap, modelConstants.theSurroundLconePlusMconePSF, 'same');

    % The composite visual RF cone map
    theVisualRF = theVisualRFcenter - theVisualRFsurround;
end